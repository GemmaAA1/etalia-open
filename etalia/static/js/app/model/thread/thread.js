define([
    'app',
    'app/model/user/user',
    'app/model/thread/state',
    'app/collection/thread/post',
    'app/collection/user/user'
], function (App) {

    App.Model.Thread = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/thread/threads',

        defaults: {
            link: null,
            type: 2,
            privacy: 1,
            title: null,
            content: null,
            state: null,
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
                type: App.Backbone.HasOne,
                key: 'paper',
                relatedModel: App.Model.Paper,
                includeInJSON: 'link'
            },
            {
                type: App.Backbone.HasOne,
                key: 'state',
                relatedModel: App.Model.State,
                includeInJSON: false,
                reverseRelation: {
                    key: 'thread',
                    type: App.Backbone.HasOne,
                    includeInJSON: 'link'
                }
            },
            {
                type: App.Backbone.HasMany,
                key: 'posts',
                relatedModel: App.Model.Post,
                collectionType: App.Collection.Posts,
                includeInJSON: false,
                reverseRelation: {
                    key: 'thread',
                    type: App.Backbone.HasOne,
                    includeInJSON: 'link'
                }
            },
            {
                type: App.Backbone.HasMany,
                key: 'members',
                relatedModel: App.Model.User,
                collectionType: App.Collection.Users,
                includeInJSON: false
            }
        ],

        getState: function() {
            var state = this.get('state');
            if (!state) {
                state = new App.Model.State({
                    user: App.Model.User.getCurrent(),
                    thread: this
                });
                this.set({state: state});
            }
            return state;
        }
    });

    App.Model.Thread.createNew = function() {
        return new App.Model.Thread({
            user: App.Model.User.getCurrent()
        });
    };

    return App.Model.Thread;
});
