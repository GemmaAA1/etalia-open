define([
    'app',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/feed/list'
], function (App) {

    /*var Router = App.Backbone.Router.extend({
        routes: {
            "":    "list",
            ":id": "detail"
        },

        list: function() {
            console.log('list');
        },

        detail: function(id) {
            console.log('detail: ' + id);
        }
    });

    var router = new Router();

    App.Backbone.history.start({
        pushState: true,
        root: '/my-feeds/'
    });

    App.$(document).on('click', '[data-router]', function(e) {
        e.preventDefault();

        router.navigate(App.$(e.currentTarget).attr('href'), true);

        return false;
    });*/

    var controlsView, filtersView, tabsView;

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
                    added: 0,
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
                    added: 0,
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

    var view = App.View.Feed.List.create({
        el: '#feed-container',
        controlsView: controlsView,
        tabsView: tabsView,
        filtersView: filtersView
    });

    App.Layout.initRightFlap();

    view.onTabsContextChange();

    App.init();
});
