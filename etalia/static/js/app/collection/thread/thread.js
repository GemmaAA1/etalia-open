define(['app', 'app/model/thread/thread'], function (App) {

    /*var defaults = {
        query: {
            view: 'nested'
        }
    };*/

    return App.Collection.Threads = App.Backbone.PageableCollection.extend({
        url: App.config.api_root + '/thread/threads',
        mode: 'infinite',
        model: App.Model.Thread,

        initialize: function(options) {

            /*console.log(options.query);
            console.log(options);

            App.defaults(options, defaults);

            console.log(options);*/

            this.queryParams = options.query;
        }
    });

});
