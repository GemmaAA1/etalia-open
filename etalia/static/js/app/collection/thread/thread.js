define(['app', 'app/model/thread/thread'], function (App) {

    var defaults = {
        query: {
            view: 'nested'
        }
    };

    return App.Collection.Threads = App.Backbone.PageableCollection.extend({
        url: App.config.api_root + '/thread/threads',
        mode: 'infinite',
        model: App.Model.Thread,

        initialize: function(options) {
            App._.defaults(options, defaults);

            this.queryParams = options.query;
        }
    });

});
