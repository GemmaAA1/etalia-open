define([
    'app',
    'text!app/templates/paper/list.hbs',
    'app/view/detail',
    'app/model/library/paper',
    'app/view/list',
    'app/view/ui/modal',
    'app/view/paper/detail',
    'app/view/paper/thumb'
], function (App, template, Detail) {

    var $window = $(window),
        $document = $(document);

    App.View.Paper = App.View.Paper || {};

    App.View.Paper.List = App.Backbone.View.extend({
        tagName: 'div',

        template: App.Handlebars.compile(template),

        altmetricTimeout: 0,

        controlsView: null,
        tabsView: null,
        filtersView: null,
        listView: null,

        events: {
            "click #paper-next-page": "onNextPageClick"
        },

        initialize: function (options) {
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
                this.listenTo(this.tabsView, "context-change", this._loadFilters);
            }

            // Filters (right flap)
            if (!options.filtersView) {
                throw 'options.filtersView is mandatory';
            }
            this.filtersView = options.filtersView;
            this.listenTo(this.filtersView, "loaded", this.load);
            this.listenTo(this.filtersView, "context-change", this.load);

            // Detail
            this.listenTo(this, "model:detail", this.openDetail);

            App._.bindAll(this,
                'onPaperPin', 'onPaperUnpin',
                'onPaperBan', 'onPaperUnban',
                'onPaperAdd', 'onPaperTrash',
                'onWindowScroll');

            App.on('etalia.paper.pin', this.onPaperPin);
            App.on('etalia.paper.unpin', this.onPaperUnpin);
            App.on('etalia.paper.ban', this.onPaperBan);
            App.on('etalia.paper.unban', this.onPaperUnban);
            App.on('etalia.paper.add', this.onPaperAdd);
            App.on('etalia.paper.trash', this.onPaperTrash);
        },

        remove: function() {
            App.off('etalia.paper.pin', this.onPaperPin);
            App.off('etalia.paper.unpin', this.onPaperUnpin);
            App.off('etalia.paper.ban', this.onPaperBan);
            App.off('etalia.paper.unban', this.onPaperUnban);
            App.off('etalia.paper.add', this.onPaperAdd);
            App.off('etalia.paper.trash', this.onPaperTrash);

            $window.off('scroll', this.onWindowScroll);

            App.Backbone.View.prototype.remove.apply(this, arguments);
        },

        render: function () {
            //App.log('PaperListView::render');

            this.$el.html(this.template({}));

            // Thumbs list
            this.listView = new App.View.List.create({}, {
                $target: this.$('div[data-list-placeholder]')
            });
            this.pushSubView(this.listView);

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

            this.collection = new App.Model.PageablePapers({
                query: App._.extend({},
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
                        that.$('#paper-next-page').show();
                        $window.on('scroll', that.onWindowScroll);
                    } else {
                        that.$('#paper-next-page').hide();
                    }
                    //App.Layout.setAvailable();

                    that._renderAltmetricBadges();
                });
        },

        onWindowScroll: function() {
            //App.log('PaperListView::onWindowScroll');
            var delta = $document.height() - $window.height() - $window.scrollTop();
            if (delta < 60) {
                this._loadNextPage();
            }
        },

        onNextPageClick: function (e) {
            e.preventDefault();

            this._loadNextPage();
        },

        _renderAltmetricBadges: function() {
            // TODO this is ugly : should be triggered once the last thumb has been rendered
            if (this.altmetricTimeout) {
                clearTimeout(this.altmetricTimeout);
            }
            this.altmetricTimeout = setTimeout(_altmetric_embed_init, 250);
        },

        _loadNextPage: function() {
            $window.off('scroll', this.onWindowScroll);

            var that = this;
            if (this.collection.hasNextPage()) {
                //App.Layout.setBusy();

                this.collection.getNextPage(null)
                    .then(function () {
                        if (that.collection.hasNextPage()) {
                            that.$('#paper-next-page').show();
                            $window.on('scroll', that.onWindowScroll);
                        } else {
                            that.$('#paper-next-page').hide();
                        }
                        //App.Layout.setAvailable();

                        that._renderAltmetricBadges();
                    });
            }
        },

        _loadFilters: function() {
            this.listView.clear();
            this.$('#paper-next-page').show();
            this.filtersView.load(
                App.config.api_root + '/library/my-papers/filters/',
                App._.extend({},
                    this.controlsView.getContext(),
                    this.tabsView.getContext()
                )
            );
        },

        onPaperPin: function(paper) {
            this._renderAltmetricBadges();

            this.tabsView.setTabCount('paper:pins', 1, true);
        },

        onPaperUnpin: function(paper) {
            var tabName = this.tabsView.getActiveTab().name;
            if (tabName == 'paper:pins') {
                this.onModelRemove(paper);
            } else {
                this._renderAltmetricBadges();
            }
            this.tabsView.setTabCount('paper:pins', -1, true);
        },

        onPaperBan: function(paper) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['feed:papers', 'paper:papers', 'paper:pins'].indexOf(tabName)) {
                this.onModelRemove(paper);
                Detail.close();
            } else {
                this._renderAltmetricBadges();
            }
            this.tabsView.setTabCount('feed:papers', -1, true);
            this.tabsView.setTabCount('paper:papers', -1, true);
        },

        onPaperUnban: function(paper) {
            /*var tabName = this.tabsView.getActiveTab().name;
            if (tabName == 'paper:trash') {
                this.onModelRemove(paper);
            } else {
                this._renderAltmetricBadges();
            }*/
            // TODO ???
        },

        onPaperAdd: function(paper) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['feed:papers', 'paper:trash'].indexOf(tabName)) {
                this.onModelRemove(paper);
            } else {
                this._renderAltmetricBadges();
            }
            this.tabsView.setTabCount('feed:papers', -1, true);
            this.tabsView.setTabCount('paper:trash', -1, true);
            this.tabsView.setTabCount('paper:papers', 1, true);
        },

        onPaperTrash: function(paper) {
            var tabName = this.tabsView.getActiveTab().name;
            if (0 <= ['paper:papers', 'paper:pins'].indexOf(tabName)) {
                this.onModelRemove(paper);
            } else {
                this._renderAltmetricBadges();
            }
            this.tabsView.setTabCount('paper:papers', -1, true);
            //this.tabsView.setTabCount('paper:pins', -1, true);
            this.tabsView.setTabCount('paper:trash', 1, true);
        },

        onCollectionAdd: function (model) {
            //App.log('PaperListView::onCollectionAdd');

            // Render the thumb
            var options = {
                id: 'paper-thumb-' + model.get('id'),
                model: model,
                list: this
            };
            if (this.tabsView.getActiveTab().actions) {
                options.buttons = this.tabsView.getActiveTab().actions;
            }
            var thumbView = new App.View.Paper.Thumb(options);

            this.listView.addThumbView(thumbView);
        },

        onCollectionRemove: function (model) {
            //App.log('PaperListView::onCollectionRemove');

            this.listView.removeThumbById('paper-thumb-' + model.get('id'));
        },

        onCollectionUpdate: function() {
            //App.log('PaperListView::onCollectionUpdate');

            this.listView.showEmptyMessage('No paper found.');
        },

        clearList: function() {
            this.listView.hideEmptyMessage().clear();
        },

        openDetail: function (model) {
            App.log('PaperListView::onModelDetail');

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
                view: new App.View.Paper.Detail(options)
            });
            detailModel.setCenterButton({
                icon: 'close',
                title: 'Back to papers list',
                callback: function() {
                    Detail.close();
                }
            });

            var prevPaper = 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null;
            if (prevPaper) {
                detailModel.setLeftButton({
                    title: 'Previous paper',
                    caption: prevPaper.get('title'),
                    callback: function() {
                        that.openDetail(prevPaper);
                    }
                });
            }

            var nextPaper = this.collection.fullCollection.at(modelIndex + 1);
            if (nextPaper) {
                detailModel.setRightButton({
                    title: 'Next paper',
                    caption: nextPaper.get('title'),
                    callback: function() {
                        that.openDetail(nextPaper);
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
            App.log('PaperListView::onModelRemove');

            this.collection.remove(model);
            this.collection.fullCollection.remove(model);
        }
    });


    App.View.Paper.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Paper.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Paper.List;
});
