define([
    'app',
    'text!app/templates/thread/list.hbs',
    'app/model/thread/thread',
    'app/view/list',
    'app/view/ui/modal',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/view/thread/invite/treat',
    'app/view/thread/form-create'
], function (App, template) {

    var $window = $(window),
        $document = $(document);

    App.View.Thread = App.View.Thread || {};

    var ListControls = App.Backbone.Model.extend({
        defaults: {
            not_published: false,
            type_paper: false,
            type_question: false,
            private: false
        },

        getContext: function() {
            var context = {};

            if (this.get('not_published')) {
                context.published = 0;
            }
            if (this.get('type_paper') || this.get('type_question')) {
                context.type = [];
                if (this.get('type_paper')) {
                    context.type.push(App.Model.Thread.TYPE_PAPER);
                }
                if (this.get('type_question')) {
                    context.type.push(App.Model.Thread.TYPE_QUESTION);
                }
            }
            if (this.get('private')) {
                context.private = 1;
            }

            return context;
        }
    });

    App.View.Thread.List = App.Backbone.View.extend({
        tagName: 'div',

        template: App.Handlebars.compile(template),

        newButton: false,
        invitesButton: false,

        listControls: null,
        controlsView: null,
        tabsView: null,
        filtersView: null,
        listView: null,

        invites: null,

        events: {
            "click #thread-create-modal": "onCreateModalClick",
            "click #thread-invites-modal": "onInvitesModalClick",

            "click #thread-toggle-not-published": "onToggleNotPublishedClick",
            "click #thread-toggle-type-paper": "onToggleTypePaperClick",
            "click #thread-toggle-type-question": "onToggleTypeQuestionClick",
            "click #thread-toggle-private": "onTogglePrivateClick",

            "click #thread-next-page": "onNextPageClick"
        },

        initialize: function (options) {
            this.newButton = options.newButton ? options.newButton : false;
            this.invitesButton = options.invitesButton ? options.invitesButton : false;

            // List controls
            this.listControls = new ListControls();
            this.listenTo(this.listControls, "change", this.onListControlsChange);

            // Controls (top bar)
            if (!options.controlsView) {
                throw 'options.controlsView is mandatory';
            }
            this.controlsView = options.controlsView;
            this.listenTo(this.controlsView, "context-change", this._loadFilters);

            // Tabs
            if (!options.tabsView) {
                throw 'options.tabsView is mandatory';
            }
            this.tabsView = options.tabsView;
            if (!options.silentTabs) {
                this.listenTo(this.tabsView, "context-change", this.onTabsContextChange);
            }

            // Filters (right flap)
            if (!options.filtersView) {
                throw 'options.filtersView is mandatory';
            }
            this.filtersView = options.filtersView;
            this.listenTo(this.filtersView, "loaded", this.load);
            this.listenTo(this.filtersView, "context-change", this.load);

            App._.bindAll(this,
                'onThreadPin', 'onThreadUnpin',
                'onThreadBan', 'onThreadUnban',
                'onThreadJoin', 'onThreadLeave',
                'onWindowScroll');

            App.on('etalia.thread.pin', this.onThreadPin);
            App.on('etalia.thread.unpin', this.onThreadUnpin);
            App.on('etalia.thread.ban', this.onThreadBan);
            App.on('etalia.thread.unban', this.onThreadUnban);
            App.on('etalia.thread.join', this.onThreadJoin);
            App.on('etalia.thread.leave', this.onThreadLeave);
        },

        remove: function() {
            App.off('etalia.thread.pin', this.onThreadPin);
            App.off('etalia.thread.unpin', this.onThreadUnpin);
            App.off('etalia.thread.ban', this.onThreadBan);
            App.off('etalia.thread.unban', this.onThreadUnban);
            App.off('etalia.thread.join', this.onThreadJoin);
            App.off('etalia.thread.leave', this.onThreadLeave);

            $window.off('scroll', this.onWindowScroll);

            App.Backbone.View.prototype.remove.apply(this, arguments);
        },

        render: function () {
            App.log('ThreadListView::render');

            this.$el.html(this.template({
                newButton: this.newButton,
                invitesButton: this.invitesButton
            }));

            // Thumbs list
            this.listView = new App.View.List.create({}, {
                $target: this.$('div[data-list-placeholder]')
            });
            this.pushSubView(this.listView);

            this._updateListControlsVisibility();
            this._updateInvitesButton();

            this._loadFilters();

            return this;
        },

        load: function() {
            //App.Layout.setBusy();

            $window.off('scroll', this.onWindowScroll);
            this.clearList();

            if (this.collection) {
                this.collection.fullCollection.off("add", this.onCollectionAdd);
                this.collection.fullCollection.off("remove", this.onCollectionRemove);
                this.collection.fullCollection.off("reset", this.onCollectionUpdate);
                this.collection.fullCollection.off("update", this.onCollectionUpdate);
                this.collection.fullCollection.reset();
                this.collection.reset();
            }

            this.collection = new App.Model.PageableThreads({
                query: App._.extend({},
                    this.listControls.getContext(),
                    this.controlsView.getContext(),
                    this.tabsView.getContext(),
                    this.filtersView.getContext()
                )
            });

            this.collection.fullCollection.on("add", this.onCollectionAdd, this);
            this.collection.fullCollection.on("remove", this.onCollectionRemove, this);
            this.collection.fullCollection.on("reset", this.onCollectionUpdate, this);
            this.collection.fullCollection.on("update", this.onCollectionUpdate, this);

            var that = this;
            this.collection.fetch()
                .then(function(data) {
                    that.tabsView.setTabCount(null, data.count);
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                        $window.on('scroll', that.onWindowScroll);
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                    //App.Layout.setAvailable();
                });
        },

        onWindowScroll: function() {
            //App.log('ThreadListView::onWindowScroll');
            var delta = $document.height() - $window.height() - $window.scrollTop();
            if (delta < 60) {
                this._loadNextPage();
            }
        },

        onNextPageClick: function (e) {
            e.preventDefault();

            this._loadNextPage();
        },

        _loadNextPage: function() {
            $window.off('scroll', this.onWindowScroll);

            var that = this;
            if (this.collection.hasNextPage()) {
                //App.Layout.setBusy();

                this.collection.getNextPage(null)
                    .then(function () {
                        if (that.collection.hasNextPage()) {
                            that.$('#thread-next-page').show();
                            $window.on('scroll', that.onWindowScroll);
                        } else {
                            that.$('#thread-next-page').hide();
                        }
                        //App.Layout.setAvailable();
                    });
            }
        },

        _loadFilters: function() {
            this.listView.clear();
            this.$('#thread-next-page').show();
            this.filtersView.load(
                App.config.api_root + '/thread/my-threads/filters/',
                App._.extend({},
                    this.listControls.getContext(),
                    this.controlsView.getContext(),
                    this.tabsView.getContext()
                )
            );
        },

        _updateListControlsVisibility: function() {
            var tab = this.tabsView.getActiveTab();

            switch(tab.name) {
                case 'thread:threads':
                    this.$('#thread-create-modal').css({display: 'inline-block'});
                    var $invites = this.$('#thread-invites-modal');
                    if (0 < $invites.data('count')) {
                        $invites.css({display: 'inline-block'});
                    }
                    break;
                case 'thread:pins':
                    this.$('#thread-create-modal, #thread-invites-modal').hide();
                    break;
                case 'thread:left':
                    this.$('#thread-create-modal, #thread-invites-modal').hide();
                    break;
            }
        },

        _updateInvitesButton: function(force) {
            var that = this,
                title = '',
                $button = this.$('#thread-invites-modal').hide();

            App.currentUserInvites
                .update(force)
                .then(function(invites) {
                    // Count
                    $button.data('count', invites.length).find('span.count').text(invites.length);

                    // Title
                    if (0 == invites.length) {
                        title = 'No pending invitation.';
                    } else if(1 == invites.length) {
                        title = '1 pending invitation';
                    } else {
                        title = invites.length + ' pending invitations';
                    }

                    // Visibility
                    $button.attr('title', title);
                    if ((0 < invites.length) && (that.tabsView.getActiveTab().name == 'thread:threads')) {
                        $button.show();
                    }
                });
        },

        onTabsContextChange: function() {
            this._updateListControlsVisibility();
            this._loadFilters();
        },

        onThreadPin: function(thread) {
            this.tabsView.setTabCount('thread:pins', 1, true);
        },

        onThreadUnpin: function(thread) {
            var tabName = this.tabsView.getActiveTab().name;
            if (tabName == 'thread:pins') {
                this.onModelRemove(thread);
            }
            this.tabsView.setTabCount('thread:pins', -1, true);
        },

        onThreadBan: function(thread) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['feed:threads', 'thread:threads', 'thread:pins'].indexOf(tabName)) {
                this.onModelRemove(thread);
            }
            this.tabsView.setTabCount('feed:threads', -1, true);
            this.tabsView.setTabCount('paper:threads', -1, true);
        },

        onThreadUnban: function(thread) {
            /*var tabName = this.tabsView.getActiveTab().name;
            if (tabName == 'left') {
                this.onModelRemove(thread);
            }*/
            // TODO ???
        },

        onThreadJoin: function(thread) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['feed:threads', 'thread:left'].indexOf(tabName)) {
                this.onModelRemove(thread);
            }
            this.tabsView.setTabCount('feed:threads', -1, true);
            this.tabsView.setTabCount('thread:left', -1, true);
            this.tabsView.setTabCount('thread:threads', 1, true);
        },

        onThreadLeave: function(thread) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['thread:threads', 'thread:pins'].indexOf(tabName)) {
                this.onModelRemove(thread);
            }
            this.tabsView.setTabCount('thread:threads', -1, true);
            //this.tabsView.setTabCount('thread:pins', -1, true);
            this.tabsView.setTabCount('thread:left', 1, true);
        },

        onCollectionAdd: function (model) {
            //App.log('ThreadListView::onCollectionAdd');

            // Render the thumb
            var options = {
                id: 'thread-thumb-' + model.get('id'),
                model: model,
                list: this
            };
            if (this.tabsView.getActiveTab().actions) {
                options.buttons = this.tabsView.getActiveTab().actions;
            }
            var thumbView = new App.View.Thread.Thumb(options);

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            //App.log('ThreadListView::onCollectionRemove');

            this.listView.removeThumbById('thread-thumb-' + model.get('id'));
        },

        onCollectionUpdate: function() {
            //App.log('ThreadListView::onCollectionUpdate');

            this.listView.showEmptyMessage('No thread found.');
        },

        clearList: function() {
            this.listView.hideEmptyMessage().clear();
        },

        onModelRemove: function (model) {
            App.log('ThreadListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        },

        onCreateModalClick: function(e) {
            e.preventDefault();

            var that = this,
                form = App.View.Thread.CreateForm.create(),
                modal = new App.View.Ui.Modal({
                    title: 'Start a new thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                form.model.save(null, {
                    wait: true,
                    validate: false,
                    success: function () {
                        that.collection.add(form.model, {at: 0});
                        modal.close();
                        //that.openDetail(form.model);
                        // TODO navigate to detail ...
                    },
                    error: function () {
                        // TODO
                        App.log('Error', arguments);
                    }
                });
            });

            modal.on('shown', function () {
                form.postRender();
            });

            form.once('cancel', function () {
                modal.close();
            });

            modal.once('hidden', function () {
                form = null;
                modal = null;
            });

            modal.render();
        },

        onInvitesModalClick: function(e) {
            e.preventDefault();

            var that = this,
                index = 0,
                invites = null,
                view = null,
                modal = null;

            function openInvite(invite) {
                var from = invite.get('from_user'),
                    title = from.get('first_name') + ' ' + from.get('last_name') + ' invites you to this thread';

                view = new App.View.Thread.InviteTreat({
                    model: invite
                });

                view.once('declined', function() {
                    index++;
                    if (index < invites.length) {
                        openInvite(invites[index]);
                    } else {
                        modal.close();
                    }
                });

                view.once('accepted', function() {
                    App.Layout.setBusy();

                    var thread = invite.get('thread');
                    thread.fetch()
                        .done(function() {
                            // Add thread to list
                            that.collection.fullCollection.add(thread, {at: 0});

                            App.Layout.setAvailable();

                            // Close modal
                            modal.close();

                            // Open detail
                            //that.openDetail(thread);
                            // TODO navigate to detail ...
                        })
                        .fail(function() {
                            App.log('Errors', 'Failed to fetch invite thread.');
                        });
                });

                if (modal) {
                    modal.updateTitle(title).updateContent(view);
                } else {
                    modal = new App.View.Ui.Modal({
                        title: title,
                        content: view,
                        footer: false
                    });

                    modal.once('hidden', function () {
                        view = null;
                        modal = null;
                        that._updateInvitesButton();
                    });

                    modal.render();
                }
            }

            App.currentUserInvites
                .get()
                .then(function(collection) {
                    invites = collection;
                    if (0 == invites.length) {
                        App.log('Errors', 'Current user has no pending invite');
                        return;
                    }
                    openInvite(invites[index]);
                });
        },

        onListControlsChange: function() {
            this._loadFilters();

            this.$('#thread-toggle-not-published')
                .toggleClass('active', this.listControls.get('not_published'));
            this.$('#thread-toggle-type-paper')
                .toggleClass('active', this.listControls.get('type_paper'));
            this.$('#thread-toggle-type-question')
                .toggleClass('active', this.listControls.get('type_question'));
            this.$('#thread-toggle-private')
                .toggleClass('active', this.listControls.get('private'));
        },

        onToggleNotPublishedClick: function(e) {
            e.preventDefault();

            this.listControls.set('not_published', !this.listControls.get('not_published'));
        },

        onToggleTypePaperClick: function(e) {
            e.preventDefault();

            var paperActive = this.listControls.get('type_paper'),
                questionActive = this.listControls.get('type_question');

            if (paperActive) {
                paperActive = false;
            } else {
                if (questionActive) {
                    questionActive = false;
                }
                paperActive = true;
            }

            this.listControls.set({
                type_question: questionActive,
                type_paper: paperActive
            });
        },

        onToggleTypeQuestionClick: function(e) {
            e.preventDefault();

            var paperActive = this.listControls.get('type_paper'),
                questionActive = this.listControls.get('type_question');

            if (questionActive) {
                questionActive = false;
            } else {
                if (paperActive) {
                    paperActive = false;
                }
                questionActive = true;
            }

            this.listControls.set({
                type_question: questionActive,
                type_paper: paperActive
            });
        },

        onTogglePrivateClick: function(e) {
            e.preventDefault();

            this.listControls.set('private', !this.listControls.get('private'));
        }
    });


    App.View.Thread.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Thread.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Thread.List;
});
