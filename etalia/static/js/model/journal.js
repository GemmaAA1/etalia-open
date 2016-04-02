define([
    'app/config',
    'backbone',
    'backbone-relational'
], function (Config, Backbone) {

    return Backbone.RelationalModel.extend({
        urlRoot: Config.api_root + '/journals',

        defaults: {
            title: null,
            short_title: null
        },

        schema: {}
    });
});
