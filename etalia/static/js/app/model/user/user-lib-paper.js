define(['app', 'app/model/library/paper'], function (App) {

    return App.Model.UserLibPaper = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/user-lib-papers',

        defaults: {
            date_created: null,
            authored: false
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'paper',
                relatedModel: App.Model.Paper,
                includeInJSON: 'link',
                autoFetch: true
            }
        ],

        schema: {}
    });
});
