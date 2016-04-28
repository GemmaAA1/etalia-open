define([
    'app',
    'app/model/user/user',
    'app/collection/thread/comment'
], function (App) {

    App.Model.Post = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/thread/posts',

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
            },
            {
                type: App.Backbone.HasMany,
                key: 'comments',
                relatedModel: App.Model.Comment,
                collectionType: App.Collection.Comments,
                includeInJSON: false,
                reverseRelation: {
                    key: 'post',
                    type: App.Backbone.HasOne,
                    includeInJSON: 'link'
                }
            }
        ]
    });

    App.Model.Post.createNew = function(thread) {
        if (!thread) {
            throw 'thread argument is mandatory';
        }
        return new App.Model.Post({
            thread: thread,
            user: App.getCurrentUser()
        });
    };

    return App.Model.Post;
});
