define(['app', 'app/model/user/user-lib', 'app/model/user/avatar'], function (App) {

    var user_path = '/user/users',
        relationship_path = '/user/relationships';

    App.Model.User = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + user_path,

        defaults: {
            email: null,
            first_name: null,
            last_name: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user_lib',
                relatedModel: App.Model.UserLib,
                includeInJSON: false, //'link',
                autoFetch: false,
                reverseRelation: {
                    key: 'user',
                    type: App.Backbone.HasOne,
                    includeInJson: 'link'
                }
            },
            {
                type: App.Backbone.HasMany,
                key: 'avatars',
                relatedModel: App.Model.Avatar,
                includeInJSON: false, //'link',
                autoFetch: false/*,
                reverseRelation: {
                    key: 'user',
                    type: App.Backbone.HasOne,
                    includeInJson: false
                }*/
            }
        ],

        schema: {
            email: {type: 'Text', validators: ['required']},
            first_name: {type: 'Text', validators: ['required']},
            last_name: {type: 'Text', validators: ['required']}
            //photo_url: {type: 'Text', validators: ['required']}
        }
    });


    App.Model.Relationship = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + relationship_path,

        defaults: {
            status: null
        },

        relations: [
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
        ]
    });

    App.Model.Relationship.STATUS_FOLLOWING = 1;
    App.Model.Relationship.STATUS_BLOCKED = 2;

    App.Model.Relationship.validateStatus = function (status) {
        return !!(status === App.Model.Relationship.STATUS_FOLLOWING
        || status === App.Model.Relationship.STATUS_BLOCKED);
    };


    App.Model.Users = App.Backbone.Collection.extend({
        url: App.config.api_root + user_path,
        model: App.Model.User
    });

    App.Model.Relationships = App.Backbone.Collection.extend({
        url: App.config.api_root + relationship_path,
        model: App.Model.Relationship
    });


    var currentUser = null;

    /**
     * Returns the current (authenticated) user.
     *
     * @return UserModel
     * @todo Promise ?
     */
    App.getCurrentUser = function () {
        if (!currentUser) {
            var id = parseInt(App.$('body').data('user-id'));
            if (!id) {
                throw 'Failed to determine user id';
            }

            currentUser = App.Model.User.find(id);
            if (!currentUser) {
                currentUser = new App.Model.User({id: id});
                currentUser.fetch()
                    .fail(function() {
                        // Redirect to home
                        window.location.href = '/';
                    });
            }

            var relationships, relationshipsXhr;

            currentUser.getRelationships = function (fetch) {
                if (!relationships) {
                    relationships = new App.Model.Relationships();
                }

                fetch = fetch || (relationshipsXhr ? false : true);

                return new Promise(function (resolve, reject) {
                    if (!relationshipsXhr || fetch) {
                        relationshipsXhr = relationships.fetch({
                            data: {'from-user': currentUser.get('id')}
                        });
                    }
                    relationshipsXhr
                        .done(function () {
                            resolve(relationships);
                        })
                        .fail(function () {
                            reject('Failed to fetch ' + type + ' relationship.');
                        });
                });
            };

            currentUser.findRelationship = function (toUser, status) {
                if (!App.Model.Relationship.validateStatus(status)) {
                    throw 'Unexpected status';
                }
                var relationShips = currentUser.getRelationships(),
                    //fromUserId = currentUser.get('id'),
                    toUserId = App.getProperty(toUser, 'id');

                return new Promise(function (resolve, reject) {
                    relationShips
                        .then(function (relationships) {
                            resolve(relationships.find(function (relationship) {
                                return relationship.get('status') === status
                                    //&& relationship.get('from_user').get('id') === fromUserId
                                    && relationship.get('to_user').get('id') === toUserId;
                            }));
                        })
                        .catch(function (error) {
                            reject(error)
                        });
                });

            };

            currentUser.isFollowed = function (user) {
                if (!user) {
                    throw 'Expected user as first argument';
                }

                var findFollowing = currentUser.findRelationship(user, App.Model.Relationship.STATUS_FOLLOWING);
                return new Promise(function (resolve, reject) {
                    findFollowing
                        .then(function (relationship) {
                            if (relationship) {
                                resolve(true);
                            } else {
                                resolve(false);
                            }
                        })
                        .catch(function (error) {
                            reject(error)
                        });
                });
            };

            currentUser.isBlocked = function (user) {
                if (!user) {
                    throw 'Expected user as first argument';
                }

                var findBlocked = currentUser.findRelationship(user, App.Model.Relationship.STATUS_BLOCKED);
                return new Promise(function (resolve, reject) {
                    findBlocked
                        .then(function (relationship) {
                            if (relationship) {
                                resolve(true);
                            } else {
                                resolve(false);
                            }
                        })
                        .catch(function (error) {
                            reject(error)
                        });
                });
            };

            currentUser.follow = function (user) {
                return new Promise(function (resolve, reject) {
                    currentUser
                        .isFollowed(user)
                        .then(function (following) {
                            if (following) {
                                resolve(true);
                            } else {
                                currentUser
                                    .isBlocked(user)
                                    .then(function (blocked) {
                                        if (blocked) {
                                            reject('User is blocked.');
                                        } else {
                                            var relationship = new App.Model.Relationship();
                                            relationship.save({
                                                from_user: currentUser,
                                                status: App.Model.Relationship.STATUS_FOLLOWING,
                                                to_user: user
                                            }, {
                                                success: function () {
                                                    relationships.push(relationship);
                                                    user.trigger('change');
                                                    resolve(true);
                                                },
                                                error: function () {
                                                    reject('Failed to create relationship');
                                                }
                                            });
                                        }
                                    })
                                    .catch(function (error) {
                                        reject(error);
                                    });
                            }
                        });

                });
            };

            currentUser.unfollow = function (user) {
                var findRelationship = currentUser.findRelationship(user, App.Model.Relationship.STATUS_FOLLOWING);
                return new Promise(function (resolve, reject) {
                    findRelationship
                        .then(function (relationship) {
                            if (relationship) {
                                relationship.destroy({
                                    success: function () {
                                        user.trigger('change');
                                        resolve(true);
                                    },
                                    error: function () {
                                        reject('Failed to remove relationship');
                                    }
                                });
                            } else {
                                reject('Following relationship not found');
                            }
                        })
                        .catch(function (error) {
                            reject(error)
                        });
                });
            };
        }
        return currentUser;
    };

    /**
     * Handlebars helpers.
     */
    App.Handlebars.registerHelper('full_name', function(user) {
        if (!user) {
            return 'Expected user as first argument';
        }
        return App.getProperty(user, 'first_name') + " " + App.getProperty(user, 'last_name');
    });

    return App.Model.User;
});
