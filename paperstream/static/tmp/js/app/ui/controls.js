define([
    'jquery',
    'app/ui/control/search',
    'app/ui/control/timespan',
    'app/ui/control/cluster',
    'app/ui/control/pinned',
    'app/ui/control/filters'
], function($, Search, Timespan, Cluster, Pinned, Filters) {

    var controls = {
        search:   new Search().init(),
        timespan: new Timespan().init(),
        cluster:  new Cluster().init(),
        pinned:   new Pinned().init(),
        filters:  new Filters().init()
    };

    controls.getStates = function() {
        return {
            time_span:    this.timespan.getValue(),
            cluster:      this.cluster.getValue(),
            pin:          this.pinned.getValue(),
            search_query: this.search.getValue(),
            filters:      this.filters.getValue()
        };
    };

    controls.setStates = function(states) {
        // Cluster
        if (states.hasOwnProperty('cluster')) {
            this.cluster.setValue(states['cluster'] || 0); // TODO remove 'or zero' (check server side)
        }
        // Time span
        if (states.hasOwnProperty('time_span')) {
            this.timespan.setValue(states['time_span']);
        }
        // Pinned
        if (states.hasOwnProperty('pin')) {
            this.pinned.setValue(states['pin'])
        }
        // Search
        if (states.hasOwnProperty('search_query')) {
            this.search.setValue(states['search_query']);
        }
        // Filters
        if (states.hasOwnProperty('filters')) {
            this.filters.render(states['filters']);
        }

        return this;
    };

    return controls;
});
