define([
    'app',
    'app/view/paper/list',
    'app/view/thread/list'
], function (App) {

    App.View.Trend = App.View.Trend || {};

    App.View.Trend.List = App.Backbone.View.extend({
        controlsView: null,
        tabsView: null,

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
                    tabsView: this.tabsView
                },
                createOptions = {
                    $target: this.$el,
                    append: true
                };
            if (tabName == 'trend:papers') {
                //this.controlsView.cluster.enable();
                //this.controlsView.timespan.enable();
                //this.controlsView.pin.enable();
                this.listView = new App.View.Paper.List.create(viewOptions, createOptions);
            } else if (tabName == 'trend:threads') {
                //this.controlsView.cluster.disable();
                //this.controlsView.timespan.disable();
                //this.controlsView.pin.disable();
                viewOptions.privacyFilter = false;
                viewOptions.publishedFilter = false;
                this.listView = new App.View.Thread.List.create(viewOptions, createOptions);
            } else {
                throw 'Unexpected tab name';
            }
        }
    });

    App.View.Trend.List.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Trend.List(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Trend.List;
});
