define([
    'app',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/thread/list',
    'app/view/thread/form-create'
], function (App) {

    var controlsView, filtersView, tabsView;

    controlsView = new App.View.Ui.Controls.create({}, {
        $target: App.$('div[data-controls-placeholder]')
    });
    controlsView.cluster.disable();
    controlsView.timespan.disable();
    controlsView.pin.disable();

    tabsView = new App.View.Ui.Tabs.create({
        tabs: [
            {
                name: 'thread:threads',
                title: 'Threads',
                icon: 'eai-comments',
                count: 0,
                data: {
                    view: 'nested',
                    joined: 1
                    //banned: "0"
                },
                actions: {
                    pin: true,
                    leave: true
                }
            },
            {
                name: 'thread:pins',
                title: 'Pins',
                icon: 'eai-pin',
                count: 0,
                data: {
                    view: 'nested',
                    pinned: 1
                },
                actions: {
                    pin: true,
                    join: true
                }
            },
            {
                name: 'thread:left',
                title: 'Left',
                icon: 'eai-leave',
                count: 0,
                data: {
                    view: 'nested',
                    left: 1
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

    App.View.Thread.List.create({
        newButton: true,
        invitesButton: true,
        controlsView: controlsView,
        tabsView: tabsView,
        filtersView: filtersView
    }, {
        $target: App.$('div[data-list-placeholder]')
    });

    App.Layout.initRightFlap();

    App.loadPopovers();
});
