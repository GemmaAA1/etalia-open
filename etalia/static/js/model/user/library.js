define([
    'app/config',
    'backbone',
    'model/user/library-paper',
    'backbone-relational'
], function (Config, Backbone, LibraryPaperModel) {

    return Backbone.RelationalModel.extend({
        urlRoot: Config.api_root + '/user/libraries',

        idAttribute: 'pk',

        defaults: {},

        schema: {},

        relations: [
            {
                type: Backbone.HasMany,
                key: 'papers',
                relatedModel: LibraryPaperModel,

                includeInJSON: false
            }
        ]
    });
});
