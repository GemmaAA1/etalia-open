define([
    'app',
    'app/model/user/user',
    'app/model/user/user-lib'
], function (App) {

    App.Model.State = App.Backbone.RelationalModel.extend({

        urlRoot: App.config.api_root + '/thread/states',

        defaults: {
            watch: null,
            participate: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user',
                relatedModel: App.Model.User,
                includeInJSON: 'link'
            }
        ],

        initialize: function() {
            this.listenTo(this, 'change', this.onChange);
        },

        onChange: function () {
            // Bubbles state change to thread change:state
            var thread = this.get('thread');
            if (thread) {
                thread.trigger('change:state', this, thread);
            }
        },

        togglePinned: function() {
            var watch = this.get('watch');

            if (watch === null) {
                watch = App.Model.State.WATCH_PINNED;
            } else if(watch === App.Model.State.WATCH_PINNED) {
                watch = null;
            } else {
                throw 'Unexpected watch state.';
            }

            this.set({watch: watch});
            // TODO this.save({watch: watch});

            return this;
        },

        toggleBanned: function() {
            var watch = this.get('watch');

            if (watch === null) {
                watch = App.Model.State.WATCH_BANNED;
            } else if (watch === App.Model.State.WATCH_BANNED) {
                watch = null;
            } else {
                throw 'Unexpected watch state.';
            }

            this.set({watch: watch});

            return this;
        },

        join: function() {
            return this._patchAction('join');
        },

        leave: function() {
            return this._patchAction('leave');
        },

        _patchAction: function(action) {
            var model = this;

            return new Promise(function(resolve, reject) {
                function doPatch() {
                    var options = {
                        type: 'PATCH',
                        dataType: 'json',
                        url: App._.result(model, 'url') + '/' + action
                    };
                    App.Backbone.ajax(options)
                        .done(function(resp) {
                            var serverAttrs = options.parse ? model.parse(resp, options) : resp;
                            if (!model.set(serverAttrs, options)) return false;
                            model.trigger('sync', model, resp, options);
                            resolve(model);
                        })
                        .fail(function() {
                            reject('Failed to patch state with action "' + action + '".');
                        });
                }

                if (model.isNew()) {
                    model.save({}, {
                        wait:true,
                        success: function() {
                            doPatch();
                        },
                        error: function() {
                            reject('Failed to create state.');
                        }
                    });
                } else {
                    doPatch();
                }
            });
        }
    });

    App.Model.State.WATCH_PINNED = 1;
    App.Model.State.WATCH_BANNED = 2;

    App.Model.State.PARTICIPATE_LEFT = 1;
    App.Model.State.PARTICIPATE_JOINED = 2;

    return App.Model.State;
});
