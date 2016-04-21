define([
    'underscore',
    'app',
    'app/model/user/user',
    'app/model/thread/state',
    'app/collection/thread/post',
    'app/collection/user/user'
], function (_, App) {

    App.Model.Thread = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/thread/threads',

        defaults: {
            link: null,
            type: null,
            privacy: null,
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

        validate: function(attrs, options) {
            //console.log(attrs, options);
            if (options && options.validate) { // Regular validation

            } else { // Form validation
                var errors = {};
                if (attrs.hasOwnProperty('privacy') && 0 > _.indexOf([App.Model.Thread.PRIVACY_PUBLIC, App.Model.Thread.PRIVACY_PRIVATE], parseInt(attrs.privacy))) {
                    errors.privacy = 'Please select a privacy.';
                }
                if (attrs.hasOwnProperty('type')) {
                    if (0 > _.indexOf([App.Model.Thread.TYPE_QUESTION, App.Model.Thread.TYPE_PAPER], parseInt(attrs.type))) {
                        errors.type = 'Please select a type.';
                    } else if (attrs.type == App.Model.Thread.TYPE_PAPER && !attrs.paper) {
                        errors.paper = 'Please select a paper.';
                    }
                }
                if (attrs.hasOwnProperty('type') && attrs.title.length < 10) {
                    errors.title = 'Please provide a title (min. 10 chars).';
                }
                if (attrs.hasOwnProperty('content') && attrs.content.length < 50) {
                    errors.title = 'Please provide a content (min. 50 chars).';
                }
                return errors;
            }
        },

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

    App.Model.Thread.PRIVACY_PUBLIC = 1;
    App.Model.Thread.PRIVACY_PRIVATE = 2;

    App.Model.Thread.TYPE_QUESTION = 1;
    App.Model.Thread.TYPE_PAPER = 2;

    App.Model.Thread.createNew = function() {
        return new App.Model.Thread({
            user: App.Model.User.getCurrent()
        });
    };

    return App.Model.Thread;
});
