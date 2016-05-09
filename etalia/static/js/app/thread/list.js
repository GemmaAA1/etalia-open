define([
    'app',
    'app/view/ui/controls',
    'app/view/ui/filters',
    'app/view/ui/tabs',
    'app/view/thread/list',
    'app/view/thread/form-create'
], function (App) {

    var controlsView, filtersView, tabsView, listView;

    controlsView = new App.View.Ui.Controls.create({}, {
        $target: App.$('div[data-controls-placeholder]')
    });
    controlsView.cluster.disable();

    listView = App.View.Thread.List.create({}, {
        $target: App.$('div[data-list-placeholder]')
    });

    tabsView = new App.View.Ui.Tabs.create({
        tabs: [
            {
                name: 'threads',
                title: 'Threads',
                data: {
                    view: 'nested',
                    //scored: "1",
                    banned: "0"
                }
            },
            {
                name: 'pins',
                title: 'Pins',
                data: {
                    view: 'nested',
                    pinned: 1
                }
            },
            {
                name: 'left',
                title: 'Left',
                data: {
                    view: 'nested',
                    left: 1
                }
            }
        ]
    }, {
        $target: this.$('div[data-tabs-placeholder]')
    });

    filtersView = App.View.Ui.Filters.create({}, {
        $target: App.$('div[data-right-flap-placeholder]')
    });

    function loadList() {
        listView.load(App._.extend(
            controlsView.getData(),
            tabsView.getData(),
            filtersView.getData()
        ));
    }

    function loadFilters() {
        App.$.get(App.config.api_root + '/thread/threads/filters/', {
            data: App._.extend(
                controlsView.getData(),
                tabsView.getData()
            ),
            dataType: 'json'
        })
        .done(function(response) {
            if (response.hasOwnProperty('users')) {
                var group = new App.Model.Ui.FilterGroup({
                        name: 'user_id',
                        label: 'Users'
                    }),
                    entries = group.get('entries');

                App._.each(response['users'], function (user) {
                    entries.add(user);
                });

                var groups = new App.Model.Ui.FilterGroups();
                groups.add(group);

                filtersView.setGroups(groups);

                loadList();
            }
        })
        .fail(function() {
            throw 'Failed to fetch "' + options.url + '" filters.';
        });
    }

    controlsView.on('controls-change', loadFilters);
    tabsView.on('tab-selection', loadFilters);
    filtersView.on('filters-change', loadList);

    loadFilters();

    App.Layout.initRightFlap();
});
