define([
    'app',
    'app/view/detail',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/paper/list'
], function (App, Detail) {

    var listView, detailView,
        controlsView, filtersView, tabsView,
        router;

    var redirectToList = function() {
        router.navigate('my-papers/', true);
    };

    var listController = function() {
        App.setHeadTitle();

        Detail.close();

        if (listView) {
            return;
        }

        controlsView = new App.View.Ui.Controls.create({}, {
            $target: App.$('div[data-controls-placeholder]')
        });
        controlsView.cluster.disable();
        controlsView.timespan.disable();
        controlsView.pin.disable();

        tabsView = new App.View.Ui.Tabs.create({
            tabs: [
                {
                    name: 'paper:papers',
                    title: 'Papers',
                    icon: 'eai-paper',
                    count: 0,
                    data: {
                        view: 'nested',
                        added: 1
                    },
                    actions: {
                        pin: true,
                        trash: true
                    }
                },
                {
                    name: 'paper:pins',
                    title: 'Pins',
                    icon: 'eai-pin',
                    count: 0,
                    data: {
                        view: 'nested',
                        pinned: 1
                    },
                    actions: {
                        pin: true,
                        add: true,
                        trash: true
                    }
                },
                {
                    name: 'paper:trash',
                    title: 'Trash',
                    icon: 'eai-library-trash',
                    count: 0,
                    data: {
                        view: 'nested',
                        trashed: 1
                    },
                    actions: {
                        add: true
                    }
                }
            ]
        }, {
            $target: App.$('div[data-tabs-placeholder]')
        });

        filtersView = App.View.Ui.Filters.create({}, {
            $target: App.$('div[data-right-flap-placeholder]')
        });

        listView = App.View.Paper.List.create({
            controlsView: controlsView,
            tabsView: tabsView,
            filtersView: filtersView
        }, {
            $target: App.$('div[data-list-placeholder]')
        });
    };

    var detailController = function (modelClass, slug) {

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
        if (modelClass === App.Model.Thread) {
            modelDetailView = new App.View.Thread.Detail(options);
        } else if (modelClass === App.Model.Paper) {
            modelDetailView = new App.View.Paper.Detail(options);
        } else {
            throw 'Unexpected model class';
        }

        /*modelDetailView.on('close', function() {
            router.navigate('my-papers/', true);
        });*/

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
                App.setHeadTitle(model.get('title'));
                Detail.setModel(detailModel);
            });

    };


    var Router = App.Backbone.Router.extend({
        routes: {
            "my-papers/": "list",
            "papers/:id/": "paperDetail",
            "threads/:id/": "threadDetail"
        },
        list: listController,
        paperDetail: function(id) {
            detailController(App.Model.Paper, id);
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

    App.Layout.initRightFlap();

    App.init();
});
