define([
    'app',
    'app/view/detail',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/feed/list'
], function (App, Detail) {

    var listView, detailView,
        controlsView, filtersView, tabsView,
        router;

    var redirectToList = function() {
        router.navigate('my-feeds/', true);
    };

    var listController = function() {
        console.log('listController()');

        Detail.close();

        if (listView) {
            return;
        }

        controlsView = new App.View.Ui.Controls.create({}, {
            $target: App.$('div[data-controls-placeholder]')
        });
        controlsView.cluster.disable();

    tabsView = new App.View.Ui.Tabs.create({
        tabs: [
            {
                name: 'feed:papers',
                title: 'Papers',
                icon: 'eai-paper',
                count: false,
                data: {
                    view: 'nested',
                    scored: 1,
                    banned: 0,
                    type: 'stream'
                },
                actions: {
                    pin: true,
                    ban: true,
                    add: true
                }
            },
            {
                name: 'feed:trend',
                title: 'Trend',
                icon: 'eai-stats',
                count: false,
                data: {
                    view: 'nested',
                    scored: 1,
                    banned: 0,
                    type: 'trend'
                },
                actions: {
                    pin: true,
                    ban: true,
                    add: true
                }
            },
            {
                name: 'feed:threads',
                title: 'Threads',
                icon: 'eai-comments',
                count: false,
                data: {
                    view: 'nested',
                    scored: 1,
                    joined: 0,
                    banned: 0
                },
                actions: {
                    pin: true,
                    ban: true,
                    join: true
                }
            }
        ]
    }, {
        $target: this.$('div[data-tabs-placeholder]')
    });

        filtersView = App.View.Ui.Filters.create({}, {
            $target: App.$('div[data-right-flap-placeholder]')
        });

        listView = App.View.Feed.List.create({
            el: '#feed-container',
            controlsView: controlsView,
            tabsView: tabsView,
            filtersView: filtersView
        });

        listView.onTabsContextChange();
    };


    var paperDetailController = function (id) {
        console.log('paperDetailController(' + id + ')');

        if (detailView) {
            Detail.close();
            detailView.remove();
            detailView = null;
        }

        // TODO parse slug
        id = parseInt(id);
        if (!id) {
            throw 'Failed to parse paper slug';
        }

        var model = App.Model.Paper.find(id);
        if (!model) {
            model = new App.Model.Paper({id: id});
        }


        var options = {
                model: model,
                listView: listView.listView // App.Paper.List.View
            };
        if (tabsView.getActiveTab().actions) {
            options.buttons = tabsView.getActiveTab().actions;
        }

        var detailModel = new App.Model.Detail({
            view: new App.View.Paper.Detail(options)
        });
        detailModel.setCenterButton({
            icon: 'close',
            title: 'Back to papers list',
            callback: function() {
                redirectToList();
            }
        });

        /*var prevPaper = 0 < modelIndex ? this.collection.fullCollection.at(modelIndex - 1) : null;
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
        }*/

        model
            .fetch({data: {view: 'nested'}})
            .fail(function() {
                redirectToList();
            })
            .done(function() {
                //console.log('detail fetch done.');
                Detail.setModel(detailModel);
            });

    };


    var Router = App.Backbone.Router.extend({
        routes: {
            "my-feeds/": "list",
            "paper/:id": "paperDetail"
        },
        list: listController,
        paperDetail: paperDetailController
    });

    router = new Router();

    App.$(document).on('click', '[data-router]', function(e) {
        e.preventDefault();

        router.navigate(App.$(e.currentTarget).attr('href'), true);

        return false;
    });

    App.Backbone.history.start({
        pushState: true
    });

    App.Layout.initRightFlap();

    App.init();

    // TODO
    // remove App.View.Thread.List.openDetail calls.
    // remove App.View.Paper.List.openDetail calls.
});
