define(['app', 'app/model/thread/post'], function (App) {

    //var defaults = {};

    return App.Collection.Posts = App.Backbone.Collection.extend({
        url: App.config.api_root + '/thread/posts',

        model: App.Model.Post,

        initialize: function(options) {
            //options = App.defaults(defaults, options);
        }
    });

});
