define([
    'jquery',
    'app',
    'app/model/user/user-lib'
], function ($, App) {

    App.Model.User = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/users',

        defaults: {
            email: null,
            first_name: null,
            last_name: null,
            photo_url: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user_lib',
                relatedModel: App.Model.UserLib,
                includeInJSON: 'link',
                autoFetch: false,
                reverseRelation: {
                    key: 'user',
                    type: App.Backbone.HasOne,
                    includeInJson: 'link'
                }
            }
        ],

        schema: {
            email: {type: 'Text', validators: ['required']},
            first_name: {type: 'Text', validators: ['required']},
            last_name: {type: 'Text', validators: ['required']}
            //photo_url: {type: 'Text', validators: ['required']}
        },

        /*constructor: function() {
            App.Backbone.RelationalModel.apply(this, arguments);

            var that = this,
                followers, followed, blocked;

            function createRelationshipPromise(type, storage) {
                return new Promise(function(resolve, reject) {
                    if (!storage) {
                        storage = new App.Collection.Users({
                            url: App.config.api_root + '/user/users/' + this.get('id') + '/' + type
                        });
                        storage.fetch()
                            .done(function() {
                                resolve(storage);
                            })
                            .fail(function() {
                                reject('Failed to fetch ' + type + ' relationship.');
                            });
                    } else {
                        resolve(storage);
                    }
                });
            }

            this.getFollowed = function() {
                return createRelationshipPromise('following', followed);
            };

            this.isFollowed = function(user) {
                if (!user) {
                    throw 'Expected user as first argument';
                }
                this.getFollowed().then(function(followed) {

                })
            };
        } */
    });

    var currentUser = null;

    /**
     * Returns the current (authenticated) user.
     *
     * @return UserModel
     */
    App.getCurrentUser = function() {
        if (!currentUser) {
            var id = parseInt($('body').data('user-id'));
            if (!id) {
                throw 'Failed to determine user id';
            }
            currentUser = App.Model.User.find(id);
            if (!currentUser) {
                currentUser = new App.Model.User({id: id});
                currentUser.fetch();
            }
        }
        return currentUser;
    };

    return App.Model.User;
});
