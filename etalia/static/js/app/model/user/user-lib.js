define(['app', 'app/model/user/user-lib-paper'], function (App) {

    return App.Model.UserLib = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/user-libs',

        defaults: {},

        relations: [
            {
                type: App.Backbone.HasMany,
                key: 'user_lib_papers',
                relatedModel: App.Model.UserLibPaper,
                includeInJSON: 'link',
                autoFetch: true,
                reverseRelation: {
                    key: 'userlib',
                    includeInJSON: 'link'
                }
            }
        ],

        schema: {}
    });
});
