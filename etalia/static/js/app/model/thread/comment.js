define(['app', 'app/model/user/user'], function (App) {

    var path = '/thread/comments';

    App.Model.Comment = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + path,

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

    App.Model.Comments = App.Backbone.Collection.extend({
        url: App.config.api_root + path,
        model: App.Model.Comment
    });

    return App.Model.Comment;
});
