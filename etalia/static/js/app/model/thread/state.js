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

            return this.save({watch: watch}, {wait: true});
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

            return this.save({watch: watch}, {wait: true});
        },

        join: function() {
            var participate = this.get('participate');

            if (participate === null || App.Model.State.PARTICIPATE_LEFT) {
                participate = App.Model.State.PARTICIPATE_JOINED;
            } else {
                throw 'Unexpected participate state.';
            }

            return this.save({participate: participate}, {wait: true});
        },

        leave: function() {
            var participate = this.get('participate');

            if (App.Model.State.PARTICIPATE_JOINED) {
                participate = App.Model.State.PARTICIPATE_LEFT;
            } else {
                throw 'Unexpected participate state.';
            }

            return this.save({participate: participate}, {wait: true});
        }
    });

    App.Model.State.WATCH_PINNED = 1;
    App.Model.State.WATCH_BANNED = 2;

    App.Model.State.PARTICIPATE_JOINED = 1;
    App.Model.State.PARTICIPATE_LEFT = 2;

    return App.Model.State;
});
