define([
    'app',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/feed/list'
], function (App) {

    var controlsView, filtersView, tabsView;

    controlsView = new App.View.Ui.Controls.create({}, {
        $target: App.$('div[data-controls-placeholder]')
    });
    controlsView.cluster.disable();

    tabsView = new App.View.Ui.Tabs.create({
        tabs: [
            {
                name: 'papers',
                title: 'Papers',
                count: 0,
                data: {
                    view: 'nested',
                    scored: 1,
                    banned: 0,
                    type: 'stream'
                },
                actions: {
                    pin: true,
                    ban: true
                }
            },
            {
                name: 'trend',
                title: 'Trend',
                count: 0,
                data: {
                    view: 'nested',
                    scored: 1,
                    banned: 0,
                    type: 'feed'
                },
                actions: {
                    pin: true
                }
            },
            {
                name: 'threads',
                title: 'Threads',
                count: 0,
                data: {
                    view: 'nested',
                    scored: 1,
                    banned: 0
                },
                actions: {
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
});
