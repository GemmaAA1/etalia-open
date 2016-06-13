define([
    'app',
    'app/model/user/user',
    'app/model/user/user-lib'
], function (App) {

    App.Model.ThreadState = App.Backbone.RelationalModel.extend({

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
                thread.trigger('change:state', this);
            }
        },

        togglePinned: function() {
            var that = this,
                watch = this.get('watch'),
                event = null;

            if (watch === null) {
                watch = App.Model.ThreadState.WATCH_PINNED;
                event = 'etalia.thread.pin';
            } else if(watch === App.Model.ThreadState.WATCH_PINNED) {
                watch = null;
                event = 'etalia.thread.unpin';
            } else {
                throw 'Unexpected watch state.';
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App.trigger(event, that.get('thread'));
                });
        },

        toggleBanned: function() {
            var that = this,
                watch = this.get('watch'),
                event = null;

            if (watch === null) {
                watch = App.Model.ThreadState.WATCH_BANNED;
                event = 'etalia.thread.ban';
            } else if (watch === App.Model.ThreadState.WATCH_BANNED) {
                watch = null;
                event = 'etalia.thread.unban';
            } else {
                throw 'Unexpected watch state.';
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App.trigger(event, that.get('thread'));
                });
        },

        join: function() {
            var that = this,
                participate = this.get('participate');

            if (participate === null || App.Model.ThreadState.PARTICIPATE_LEFT) {
                participate = App.Model.ThreadState.PARTICIPATE_JOINED;
            } else {
                throw 'Unexpected participate state.';
            }

            return this
                .save({participate: participate}, {patch: true, wait: true})
                .done(function() {
                    App.trigger('etalia.thread.join', that.get('thread'));
                });
        },

        leave: function() {
            var that = this,
                participate = this.get('participate');

            if (App.Model.ThreadState.PARTICIPATE_JOINED) {
                participate = App.Model.ThreadState.PARTICIPATE_LEFT;
            } else {
                throw 'Unexpected participate state.';
            }

            return this
                .save({participate: participate}, {patch: true, wait: true})
                .done(function() {
                    App.trigger('etalia.thread.leave', that.get('thread'));
                });
        }
    });

    App.Model.ThreadState.WATCH_PINNED = 1;
    App.Model.ThreadState.WATCH_BANNED = 2;

    App.Model.ThreadState.PARTICIPATE_JOINED = 1;
    App.Model.ThreadState.PARTICIPATE_LEFT = 2;

    return App.Model.ThreadState;
});
