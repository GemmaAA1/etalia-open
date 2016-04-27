define(['app', 'app/model/user/user-lib'], function (App) {

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
                includeInJSON: false, //'link',
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
        }
    });


    App.Model.Relationship = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/relationships',

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


    var currentUser = null;

    /**
     * Returns the current (authenticated) user.
     *
     * @return UserModel
     * @todo to Promise
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
                currentUser.fetch();
            }

            var relationships, relationshipsXhr;

            currentUser.getRelationships = function (fetch) {
                if (!relationships) {
                    relationships = new App.Collection.Relationships();
                }

                fetch = fetch || (relationshipsXhr ? false : true);

                return new Promise(function (resolve, reject) {
                    if (!relationshipsXhr || fetch) {
                        relationshipsXhr = relationships.fetch({
                            //url: App.config.api_root + '/user/users/' + that.get('id') + '/' + type
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
                    fromUserId = currentUser.get('id'),
                    toUserId = toUser.hasOwnProperty('get') ? toUser.get('id') : toUser.id;

                return new Promise(function (resolve, reject) {
                    relationShips
                        .then(function (relationships) {
                            resolve(relationships.find(function (relationship) {
                                return relationship.get('status') === status
                                        // TODO Remove as soon as we fetch only "this" user relationships
                                    && relationship.get('from_user').get('id') === fromUserId
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

    return App.Model.User;
});
