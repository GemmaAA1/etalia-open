define([
    'app/config',
    'backbone',
    'backbone-relational'
], function (Config, Backbone) {

    return Backbone.RelationalModel.extend({
        //urlRoot: Config.api_root + '/authors',

        defaults: {
            first_name: null,
            last_name: null
        },

        schema: {}
    });
});
