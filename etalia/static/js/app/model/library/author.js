define(['app'], function (App) {

    return App.Model.Author = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/library/authors',

        defaults: {
            first_name: null,
            last_name: null
        },

        schema: {}
    });
});
