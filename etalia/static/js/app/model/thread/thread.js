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
            type: 1,
            privacy: 1,
            title: null,
            content: null,
            created: null,
            modified: null,
            published_at: null
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
                if (attrs.hasOwnProperty('title') && attrs.title.length < 10) {
                    errors.title = 'Please provide a title (min. 10 chars).';
                }
                if (attrs.hasOwnProperty('content') && attrs.content.length < 50) {
                    errors.content = 'Please provide a content (min. 50 chars).';
                }
                return errors;
            }
        },

        publish: function() {
            var options = {
                type: 'PATCH',
                dataType: 'json',
                url: _.result(this, 'url') + '/publish'
            };

            var model = this;
            App.Backbone.ajax(options).then(function(resp) {
                var serverAttrs = options.parse ? model.parse(resp, options) : resp;
                if (!model.set(serverAttrs, options)) return false;
                model.trigger('sync', model, resp, options);
            });

            return this;
        },

        getState: function() {
            var state = this.get('state');
            if (!state) {
                state = new App.Model.State({
                    user: App.getCurrentUser(),
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
            user: App.getCurrentUser()
        });
    };

    return App.Model.Thread;
});
