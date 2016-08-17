define(['app/app'], function (App) {

    var path = '/library/authors';

    App.Model.Author = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + path,

        defaults: {
            first_name: null,
            last_name: null
        },

        schema: {}
    });

    App.Model.Authors = App.Backbone.Collection.extend({
        url: App.config.api_root + path,
        model: App.Model.Author
    });

    return App.Model.Author
});
