define([
    'app',
    'text!app/templates/thread/list.hbs',
    'app/model/thread/thread',
    'app/view/list',
    'app/view/detail',
    'app/view/ui/modal',
    'app/view/thread/detail',
    'app/view/thread/thumb',
    'app/view/thread/invite/treat'
], function (App, template) {

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
            // TODO type (paper / question)
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
        detailView: null,

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
            App.Layout.setBusy();

            this.clearList();

            if (this.collection) {
                this.collection.fullCollection.off("add", this.onCollectionAdd);
                this.collection.fullCollection.off("remove", this.onCollectionRemove);
                this.collection.fullCollection.reset();
                this.collection.reset();
            }

            this.collection = new App.Model.Threads({
                query: App._.extend(
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
                .then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                    App.Layout.setAvailable();
                });
        },

        _loadFilters: function() {
            this.filtersView.load(
                App.config.api_root + '/thread/threads/filters/',
                App._.extend(
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
            var thumbView = new App.View.Thread.Thumb({
                id: 'thread-thumb-' + model.get('id'),
                model: model,
                list: this
            });

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
                detailModel = new App.Model.Detail({
                //list: this,
                //prev: 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null,
                //next: this.collection.fullCollection.at(modelIndex + 1),
                view: new App.View.Thread.Detail({
                    model: model
                })
            });
            detailModel.setCenterButton({
                icon: 'close',
                title: 'Back to threads list',
                callback: function() {
                    if (that.detailView) that.detailView.close();
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

            if (this.detailView) {
                this.detailView.clearSubViews();
                this.detailView.model.destroy();
                this.detailView.model = detailModel;
            } else {
                this.detailView = new App.View.Detail({
                    listView: this,
                    model: detailModel
                });
                this.listenTo(this.detailView, "detail:prev", this.openDetail);
                this.listenTo(this.detailView, "detail:next", this.openDetail);
            }

            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    that.detailView.render();
                });
        },

        onModelRemove: function (model) {
            App.log('ThreadListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        },

        onNextPageClick: function (e) {
            e.preventDefault();

            var that = this;
            if (this.collection.hasNextPage()) {
                this.collection.getNextPage(null).then(function() {
                    if (that.collection.hasNextPage()) {
                        that.$('#thread-next-page').show();
                    } else {
                        that.$('#thread-next-page').hide();
                    }
                });
            }
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
