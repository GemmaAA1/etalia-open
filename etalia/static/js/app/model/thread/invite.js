define([
    'app',
    'app/model/user/user',
    'app/model/thread/thread'
], function (App) {

    App.Model.Invite = App.Backbone.RelationalModel.extend({

        urlRoot: App.config.api_root + '/thread/invites',

        defaults: {
            status: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'thread',
                relatedModel: App.Model.Thread,
                includeInJSON: 'link'
            },
            {
                type: App.Backbone.HasOne,
                key: 'from_user',
                relatedModel: App.Model.User,
                includeInJSON: 'link'
            },
            {
                type: App.Backbone.HasOne,
                key: 'to_user',
                relatedModel: App.Model.User,
                includeInJSON: 'link'
            }
        ],

        validate: function(attrs, options) {
            var errors = {};

            if (options && options.validate) {
                // TODO use model validators
                if (attrs.hasOwnProperty('status') && !App.Model.Invite.isValidStatus(attrs.status)) {
                    errors.status = 'Please select a status.';
                }
                if (attrs.hasOwnProperty('thread') && !App.Validator.isRelation(attrs.thread)) {
                    errors.thread = 'Please select a thread.';
                }
                if (attrs.hasOwnProperty('from_user') && !App.Validator.isRelation(attrs.from_user)) {
                    errors.from_user = 'Please select a user.';
                }
                if (attrs.hasOwnProperty('to_user') && !App.Validator.isRelation(attrs.to_user)) {
                    errors.to_user = 'Please select a user.';
                }
            }

            if (!App._.isEmpty(errors)) {
                App.log('Errors', errors);
                return errors;
            }
        },

        accept: function() {
            var status = this.get('status');

            if (status === null || App.Model.Invite.STATUS_PENDING) {
                status = App.Model.Invite.STATUS_ACCEPTED;
            } else {
                throw 'Unexpected invite status.';
            }

            return this.save({status: status}, {patch: true, wait: true});
        },

        decline: function() {
            var status = this.get('status');

            if (status === null || App.Model.Invite.STATUS_PENDING) {
                status = App.Model.Invite.STATUS_DECLINED;
            } else {
                throw 'Unexpected invite status.';
            }

            return this.save({status: status}, {patch: true, wait: true});
        }
    });

    App.Model.Invite.validators = {
        status: function (status) {
            if(!App.Model.Invite.isValidStatus(status)) {
                return {
                    type: 'status',
                    message: 'Invalid status.'
                }
            }
        },
        thread: function (thread) {
            if(!App.Validator.isRelation(thread)) {
                return {
                    type: 'from_user',
                    message: 'Please select the thread.'
                }
            }
        },
        from_user: function (user) {
            if(!App.Validator.isRelation(user)) {
                return {
                    type: 'from_user',
                    message: 'Please select the source user.'
                }
            }
        },
        to_user: function (user) {
            if(!App.Validator.isRelation(user)) {
                return {
                    type: 'to_user',
                    message: 'Please select the target user.'
                }
            }
        }
    };

    App.Model.Invite.STATUS_PENDING  = 1;
    App.Model.Invite.STATUS_ACCEPTED = 2;
    App.Model.Invite.STATUS_DECLINED = 3;

    App.Model.Invite.isValidStatus = function (status) {
        return 0 <= App._.indexOf([
            App.Model.Invite.STATUS_PENDING,
            App.Model.Invite.STATUS_ACCEPTED,
            App.Model.Invite.STATUS_DECLINED
        ], parseInt(status));
    };

    App.Model.Invite.createNew = function(data) {
        data = App._.defaults(data, {
            status: App.Model.Invite.STATUS_PENDING,
            from_user: App.getCurrentUser()
        });

        return new App.Model.Invite(data);
    };

    App.Model.Invites = App.Backbone.Collection.extend({
        url: App.config.api_root + '/thread/invites',
        model: App.Model.Invite
    });

    var currentUserInvites = function() {
        var invites = [];

        this.update = function(force) {
            force = force || false;

            return new Promise(function(resolve, reject) {
                if (!force && 0 < invites.length) {
                    resolve(invites);
                } else {
                    var collection = new App.Model.Invites();
                    collection
                        .fetch()
                        .done(function () {
                            invites = collection.where({
                                to_user: App.getCurrentUser(),
                                status: App.Model.Invite.STATUS_PENDING
                            });
                            resolve(invites);
                        })
                        .fail(function () {
                            reject('Failed to fetch current user pending invites.');
                        });
                }
            });
        };

        this.get = function() {
            if (0 < invites.length) {
                return new Promise(function(resolve, reject) {
                    resolve(invites);
                });
            }
            return this.update();
        };
    };


    App.currentUserInvites = new currentUserInvites();

    return App.Model.Invite;
});
