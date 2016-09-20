define([
    'app',
    'app/view/detail',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/trend/list'
], function (App, Detail) {

    var listView, detailView,
        controlsView, tabsView, defaultActiveTab = 'trend:papers',
        router;

    var redirectToList = function() {
        router.navigate('trends/', true);
    };

    var listController = function() {

        //console.log('list controller');

        Detail.close();

        if (listView) {
            return;
        }

        controlsView = new App.View.Ui.Controls.create({}, {
            $target: App.$('div[data-controls-placeholder]')
        });
        controlsView.cluster.disable();
        controlsView.pin.disable();

        tabsView = new App.View.Ui.Tabs.create({
            tabs: [
                {
                    name: 'trend:papers',
                    title: 'Papers',
                    icon: 'eai-paper',
                    count: false,
                    data: {
                        view: 'nested',
                        //type: 'trend',
                        ordering: '-altmetric__score'
                    },
                    actions: {}
                },
                {
                    name: 'trend:threads',
                    title: 'Threads',
                    icon: 'eai-comments',
                    count: false,
                    data: {
                        view: 'nested'
                    },
                    actions: {}
                }
            ]
        }, {
            $target: App.$('div[data-tabs-placeholder]')
        });
        tabsView.setActiveTab(defaultActiveTab);

        /*filtersView = App.View.Ui.Filters.create({}, {
            $target: App.$('div[data-right-flap-placeholder]')
        });*/

        listView = App.View.Trend.List.create({
            el: '#trend-container',
            controlsView: controlsView,
            tabsView: tabsView
            //filtersView: filtersView
        });

        listView.onTabsContextChange();
    };

    var detailController = function (modelClass, slug) {

        //console.log('detail controller');

        if (detailView) {
            Detail.close();
            detailView.remove();
            detailView = null;
        }

        try {
            var id = parseInt(slug.match(/([a-z0-9-]+)_(\d+)/)[2]);
        } catch (e) {
            throw 'Failed to parse slug';
        }
        if (!id) {
            throw 'Failed to parse slug';
        }

        var model = modelClass.find(id);
        if (!model) {
            model = new modelClass({id: id});
        }

        var options = {model: model};
        if (tabsView && tabsView.getActiveTab().actions) {
            options.buttons = tabsView.getActiveTab().actions;
        }

        var modelDetailView;
        if (modelClass == App.Model.Thread) {
            if (!listView) {
                defaultActiveTab = 'trend:threads';
            }
            modelDetailView = new App.View.Thread.Detail(options);
        } else if (modelClass == App.Model.Paper) {
            if (!listView) {
                defaultActiveTab = 'trend:papers';
            }
            modelDetailView = new App.View.Paper.Detail(options);
        } else {
            throw 'Unexpected model class';
        }
        var detailModel = new App.Model.Detail({
            view: modelDetailView
        });
        detailModel.setCenterButton({
            icon: 'close',
            title: 'Back to papers list',
            callback: function() {
                redirectToList();
            }
        });

        model
            .fetch({data: {view: 'nested'}, silent: true})
            .fail(function() {
                redirectToList();
            })
            .done(function() {
                Detail.setModel(detailModel);
            });

    };


    var Router = App.Backbone.Router.extend({
        routes: {
            "trends/": "list",
            "papers/:slug/": "paperDetail",
            "threads/:id/": "threadDetail"
        },
        list: listController,
        paperDetail: function(slug) {
            detailController(App.Model.Paper, slug);
        },
        threadDetail: function(id) {
            detailController(App.Model.Thread, id);
        }
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

    //App.Layout.initRightFlap();

    App.init();

});
