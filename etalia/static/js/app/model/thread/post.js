define(['app', 'app/model/thread/comment'], function (App) {

    var path = '/thread/posts';

    App.Model.Post = App.Backbone.RelationalModel.extend({
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
            },
            {
                type: App.Backbone.HasMany,
                key: 'comments',
                relatedModel: App.Model.Comment,
                collectionType: App.Model.Comments,
                includeInJSON: false,
                parse: true,
                reverseRelation: {
                    key: 'post',
                    type: App.Backbone.HasOne,
                    includeInJSON: 'link'
                }
            }
        ],

        parse: function(response, options) {
            if (response.hasOwnProperty('thread') && App._.isObject(response.thread)) {
                response.thread = response.thread.link;
            }
            return response;
        },

        isOwner: function (user) {
            return this.get('user').get('id') === user.get('id');
        }
    });

    App.Model.Posts = App.Backbone.Collection.extend({
        url: App.config.api_root + path,
        model: App.Model.Post
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
