define([
    'app',
    'app/model/user/user',
    'app/model/user/user-lib',
    'app/model/library/paper'
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
            // TODO this.save({watch: watch});

            return this;
        }
    });

    App.Model.State.WATCH_PINNED = 1;
    App.Model.State.WATCH_BANNED = 2;

    App.Model.State.PARTICIPATE_LEFT = 1;
    App.Model.State.PARTICIPATE_JOINED = 2;

    return App.Model.State;
});
