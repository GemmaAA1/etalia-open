define([
    'app',
    'app/view/paper/list',
    'app/view/thread/list'
], function (App) {

    App.View.Feed = App.View.Feed || {};

    App.View.Feed.List = App.Backbone.View.extend({
        controlsView: null,
        tabsView: null,
        filtersView: null,

        listView: null,

        initialize: function (options) {
            // Tabs
            if (!options.tabsView) {
                throw 'options.tabsView is mandatory';
            }
            this.tabsView = options.tabsView;
            this.listenTo(this.tabsView, "context-change", this.onTabsContextChange);

            // Controls (top bar)
            if (!options.controlsView) {
                throw 'options.controlsView is mandatory';
            }
            this.controlsView = options.controlsView;

            // Filters (right flap)
            if (!options.filtersView) {
                throw 'options.filtersView is mandatory';
            }
            this.filtersView = options.filtersView;
        },

        onTabsContextChange: function() {
            if (this.listView) {
                this.listView.remove();
                this.listView = null;
                this.$el.empty();
            }
            var tabName = this.tabsView.getActiveTab().name,
                viewOptions = {
                    silentTabs: true,
                    controlsView: this.controlsView,
                    tabsView: this.tabsView,
                    filtersView: this.filtersView
                },
                createOptions = {
                    $target: this.$el,
                    append: true
                };
            if (0 <= ['feed:papers', 'feed:trend'].indexOf(tabName)) {
                //this.controlsView.cluster.enable();
                //this.controlsView.timespan.enable();
                //this.controlsView.pin.enable();

                this.listView = new App.View.Paper.List.create(viewOptions, createOptions);
            } else if (tabName == 'feed:threads') {
                //this.controlsView.cluster.disable();
                //this.controlsView.timespan.disable();
                //this.controlsView.pin.disable();

                this.listView = new App.View.Thread.List.create(viewOptions, createOptions);
            } else {
                throw 'Unexpected tab name';
            }
        }
    });

    App.View.Feed.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Feed.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Feed.List;
});
