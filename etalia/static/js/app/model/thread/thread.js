define([
    'app',
    'app/model/thread/state',
    'app/model/library/paper',
    'app/collection/thread/post',
    'app/collection/user/user'
], function (App) {

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

        /*constraints: {
            privacy: [
                {
                    type: 'Choice',
                    groups: ['default', 'create', 'edit'],
                    choices: [1, 2], // TODO use constants
                    message: 'Invalid privacy choice.'
                }
            ]
        },*/

        initialize: function() {
            this.listenTo(this, 'change:posts', function() {
                if (0 < this.get('posts').length) {
                    this.set('posts_count', this.get('posts').length);
                }
            });
            this.listenTo(this, 'change:members', function() {
                if (0 < this.get('members').length) {
                    this.set('members_count', this.get('members').length);
                }
            });
        },

        validate: function(attrs, options) {
            var errors = {};

            if (options && options.validate) {
                // TODO use model validators
                if (attrs.hasOwnProperty('privacy') && !App.Model.Thread.isValidPrivacy(attrs.privacy)) {
                    errors.privacy = 'Please select a privacy.';
                }
                if (attrs.hasOwnProperty('type')) {
                    if (!App.Model.Thread.isValidType(attrs.type)) {
                        errors.type = 'Please select a type.';
                    } else if (attrs.type == App.Model.Thread.TYPE_PAPER && !attrs.paper) {
                        errors.paper = 'Please select a paper.';
                    }
                }
                if (attrs.hasOwnProperty('title') && !(App._.isString(attrs.title) && 10 <= attrs.title.length)) {
                    errors.title = 'Please provide a title (min. 10 chars).';
                }
                if (attrs.hasOwnProperty('content') && !(App._.isString(attrs.content) && 50 <= attrs.content.length)) {
                    errors.content = 'Please provide a content (min. 50 chars).';
                }
            }

            if (!App._.isEmpty(errors)) {
                App.log('Errors', errors);
                return errors;
            }
        },

        publish: function() {
            var options = {
                type: 'PATCH',
                dataType: 'json',
                url: App._.result(this, 'url') + '/publish'
            };

            var model = this;
            App.Backbone.ajax(options)
                .done(function(resp) {
                    var serverAttrs = options.parse ? model.parse(resp, options) : resp;
                    if (!model.set(serverAttrs, options)) return false;
                    model.trigger('sync', model, resp, options);
                })
                .error(function() {
                    App.log('Errors', 'Failed to publish thread.');
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
        },

        isPrivate: function() {
            return this.get('privacy') === App.Model.Thread.PRIVACY_PRIVATE;
        },

        isPublic: function() {
            return this.get('privacy') === App.Model.Thread.PRIVACY_PUBLIC;
        },

        isOwner: function (user) {
            return this.get('user').get('id') === user.get('id');
        },

        isMember: function (user) {
            return this.get('members').some(function (member) {
                return member.get('id') === user.get('id');
            })
        },

        getMembersCount: function() {
            return this.getRelation('members').getCount()
        },

        getPostsCount: function() {
            return this.getRelation('posts').getCount()
        }
    });

    App.Model.Thread.validators = {
        privacy: function (privacy) {
            if(!App.Model.Thread.isValidPrivacy(privacy)) {
                return {
                    type: 'privacy',
                    message: 'Please select a privacy.'
                }
            }
        },
        type: function (type) {
            if(!App.Model.Thread.isValidType(type)) {
                return {
                    type: 'type',
                    message: 'Please select a type.'
                }
            }
        },
        paper: function (paper, formValues) {
            if (formValues['type'] == App.Model.Thread.TYPE_PAPER && !parseInt(paper)) {
                return {
                    type: 'paper',
                    message: 'Please select a paper'
                }
            }
        },
        title: function(title) {
            if (10 > String(title).length) {
                return {
                    type: 'title',
                    message: 'Title should be at least 10 characters long'
                }
            }
        }
    };

    App.Model.Thread.PRIVACY_PUBLIC = 1;
    App.Model.Thread.PRIVACY_PRIVATE = 2;

    App.Model.Thread.TYPE_QUESTION = 1;
    App.Model.Thread.TYPE_PAPER = 2;

    App.Model.Thread.isValidPrivacy = function (privacy) {
        return 0 <= App._.indexOf([
            App.Model.Thread.PRIVACY_PUBLIC,
            App.Model.Thread.PRIVACY_PRIVATE
        ], parseInt(privacy));
    };

    App.Model.Thread.isValidType = function (type) {
        return 0 <= App._.indexOf([
            App.Model.Thread.TYPE_QUESTION,
            App.Model.Thread.TYPE_PAPER
        ], parseInt(type));
    };

    App.Model.Thread.createNew = function() {
        return new App.Model.Thread({
            user: App.getCurrentUser()
        });
    };

    /**
     * Handlebars helpers
     */
    App.Handlebars.registerHelper('thread_pin_class', function() {
        if (this.state && this.state.get('watch') === App.Model.State.WATCH_PINNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('thread_ban_class', function() {
        if (this.state && this.state.get('watch') === App.Model.State.WATCH_BANNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('thread_privacy_icon', function() {
        if (this.privacy && this.privacy === App.Model.Thread.PRIVACY_PRIVATE) {
            return new App.Handlebars.SafeString('<span class="eai eai-locked"></span>');
        }
        return '';
    });
    App.Handlebars.registerHelper('thread_type_icon', function() {
        var icon = '<span class="eai eai-question"></span>';
        if (this.type && this.type === App.Model.Thread.TYPE_PAPER) {
            icon = '<span class="eai eai-paper"></span>';
        }
        return new App.Handlebars.SafeString(icon);
    });
    App.Handlebars.registerHelper('thread_type_title', function() {
        if (this.type) {
            if (this.type === App.Model.Thread.TYPE_PAPER) {
                return 'About a paper';
            } else if (this.type === App.Model.Thread.TYPE_QUESTION) {
                return 'Question';
            }
        }
        return 'Unknown';
    });

    return App.Model.Thread;
});
