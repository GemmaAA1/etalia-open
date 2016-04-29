define([
    'app',
    'app/model/user/user'
], function (App) {

    App.Model.Comment = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/thread/comments',

        defaults: {
            position: null,
            content: null,
            created: null,
            modified: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user',
                relatedModel: App.Model.User,
                includeInJSON: 'link'
            }
        ],

        isOwner: function (user) {
            return this.get('user').get('id') === user.get('id');
        }
    });

    App.Model.Comment.createNew = function() {
        return new App.Model.Comment({
            'author': App.getCurrentUser()
        });
    };

    return App.Model.Comment;
});
