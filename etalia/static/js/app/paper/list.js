define([
    'app',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/paper/list'
], function (App) {

    var controlsView, filtersView, tabsView;

    controlsView = new App.View.Ui.Controls.create({}, {
        $target: App.$('div[data-controls-placeholder]')
    });
    controlsView.cluster.disable();

    tabsView = new App.View.Ui.Tabs.create({
        tabs: [
            {
                name: 'paper:papers',
                title: 'Papers',
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
        $target: this.$('div[data-tabs-placeholder]')
    });

    filtersView = App.View.Ui.Filters.create({}, {
        $target: App.$('div[data-right-flap-placeholder]')
    });

    App.View.Paper.List.create({
        controlsView: controlsView,
        tabsView: tabsView,
        filtersView: filtersView
    }, {
        $target: App.$('div[data-list-placeholder]')
    });

    App.Layout.initRightFlap();
});
