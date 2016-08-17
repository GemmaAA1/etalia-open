define(['app/app'], function (App) {

    return App.Model.Journal = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/library/journals',

        defaults: {
            title: null,
            short_title: null,
            url: null
        },

        schema: {}
    });
});
