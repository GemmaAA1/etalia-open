define([
    'app',
    'text!app/templates/thread/list.hbs',
    'app/view/detail',
    'app/model/thread/thread',
    'app/view/list',
    'app/view/ui/modal',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/view/thread/invite/treat'
], function (App, template, Detail) {

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
            // List controls
            this.listControls = new ListControls();
            this.listenTo(this.listControls, "change", this._loadFilters);

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
            this.listenTo(this.tabsView, "context-change", this.onTabsContextChange);

            // Filters (right flap)
            if (!options.filtersView) {
                throw 'options.filtersView is mandatory';
            }
            this.filtersView = options.filtersView;
            this.listenTo(this.filtersView, "loaded", this.load);
            this.listenTo(this.filtersView, "context-change", this.load);

            // Collection
            this.listenTo(this, "model:remove", this.onModelRemove);
            this.listenTo(this, "model:detail", this.openDetail);

            this.onWindowScroll = App._.bind(this.onWindowScroll, this);
        },

        remove: function() {
            $window.off('scroll', this.onWindowScroll);
            App.Backbone.View.prototype.remove.apply(this, arguments);
        },

        render: function () {
            App.log('ThreadListView::render');

            this.$el.html(this.template({}));

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
            App.log('ThreadListView::onWindowScroll');
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
            this.filtersView.load(
                App.config.api_root + '/thread/threads/filters/',
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
                case 'threads':
                    this.$('#thread-create-modal').css({display: 'inline-block'});
                    var $invites = this.$('#thread-invites-modal');
                    if (0 < $invites.data('count')) {
                        $invites.css({display: 'inline-block'});
                    }
                    break;
                case 'pins':
                    this.$('#thread-create-modal, #thread-invites-modal').hide();
                    break;
                case 'left':
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
                    if ((0 < invites.length) && (that.tabsView.getActiveTab().name == 'threads')) {
                        $button.show();
                    }
                });
        },

        onTabsContextChange: function() {
            this._updateListControlsVisibility();
            this._loadFilters();
        },

        onCollectionAdd: function (model) {
            App.log('ThreadListView::onCollectionAdd');

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
            App.log('ThreadListView::onCollectionRemove');

            this.listView.removeThumbById('thread-thumb-' + model.get('id'));
        },

        clearList: function() {
            this.listView.clear();
        },

        openDetail: function (model) {
            App.log('ThreadListView::onModelDetail');

            var modelIndex = this.collection.fullCollection.indexOf(model);
            if (0 > modelIndex) {
                throw 'Unexpected index';
            }

            var that = this,
                options = {
                    model: model,
                    listView: this
                };
            if (this.tabsView.getActiveTab().actions) {
                options.buttons = this.tabsView.getActiveTab().actions;
            }

            var detailModel = new App.Model.Detail({
                view: new App.View.Thread.Detail(options)
            });
            detailModel.setCenterButton({
                icon: 'close',
                title: 'Back to threads list',
                callback: function() {
                    Detail.close();
                }
            });

            var prevThread = 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null;
            if (prevThread) {
                detailModel.setLeftButton({
                    title: 'Previous thread',
                    caption: prevThread.get('title'),
                    callback: function() {
                        that.openDetail(prevThread);
                    }
                });
            }

            var nextThread = this.collection.fullCollection.at(modelIndex + 1);
            if (nextThread) {
                detailModel.setRightButton({
                    title: 'Next thread',
                    caption: nextThread.get('title'),
                    callback: function() {
                        that.openDetail(nextThread);
                    }
                });
            }

            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    Detail.setModel(detailModel);
                });
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
                        that.openDetail(form.model);
                    },
                    error: function () {
                        // TODO
                        App.log('Error', arguments);
                    }
                });
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
                            that.openDetail(thread);
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

        onToggleNotPublishedClick: function(e) {
            e.preventDefault();

            var active = !this.listControls.get('not_published');
            this.listControls.set('not_published', active);

            if (active) {
                this.$('#thread-toggle-not-published').addClass('active');
            } else {
                this.$('#thread-toggle-not-published').removeClass('active');
            }
        },

        onToggleTypePaperClick: function(e) {
            e.preventDefault();

            var active = !this.listControls.get('type_paper');
            this.listControls.set('type_paper', active);

            if (active) {
                this.$('#thread-toggle-type-paper').addClass('active');
            } else {
                this.$('#thread-toggle-type-paper').removeClass('active');
            }
        },

        onToggleTypeQuestionClick: function(e) {
            e.preventDefault();

            var active = !this.listControls.get('type_question');
            this.listControls.set('type_question', active);

            if (active) {
                this.$('#thread-toggle-type-question').addClass('active');
            } else {
                this.$('#thread-toggle-type-question').removeClass('active');
            }
        },

        onTogglePrivateClick: function(e) {
            e.preventDefault();

            var active = !this.listControls.get('private');
            this.listControls.set('private', active);

            if (active) {
                this.$('#thread-toggle-private').addClass('active');
            } else {
                this.$('#thread-toggle-private').removeClass('active');
            }
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
