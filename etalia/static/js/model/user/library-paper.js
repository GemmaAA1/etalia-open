define([
    'app/config',
    'backbone',
    'model/paper',
    'backbone-relational'
], function (Config, Backbone, PaperModel) {

    return Backbone.RelationalModel.extend({
        //urlRoot: Config.api_root + '/user/libraries',

        defaults: {
            date_created: null,
            authored: false
        },

        schema: {},

        relations: [
            {
                type: Backbone.HasOne,
                key: 'paper',
                relatedModel: PaperModel,
                includeInJSON: false
            }
        ]
    });
});
