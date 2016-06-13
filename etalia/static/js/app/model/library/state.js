define([
    'app',
    'app/model/user/user',
    'app/model/user/user-lib'
], function (App) {

    App.Model.PaperState = App.Backbone.RelationalModel.extend({

        urlRoot: App.config.api_root + '/library/states',

        defaults: {
            watch: null,
            store: null
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
            // Bubbles state change to paper change:state
            var paper = this.get('paper');
            if (paper) {
                paper.trigger('change:state', this);
            }
        },

        togglePinned: function() {
            var that = this,
                watch = this.get('watch'),
                event = null;

            if (watch === null) {
                watch = App.Model.PaperState.WATCH_PINNED;
                event = 'etalia.paper.pin';
            } else if(watch === App.Model.PaperState.WATCH_PINNED) {
                watch = null;
                event = 'etalia.paper.unpin';
            } else {
                throw 'Unexpected watch state.';
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App.trigger(event, that.get('paper'));
                });
        },

        toggleBanned: function() {
            var that = this,
                watch = this.get('watch'),
                event = null;

            if (watch === null) {
                watch = App.Model.PaperState.WATCH_BANNED;
                event = 'etalia.paper.ban';
            } else if (watch === App.Model.PaperState.WATCH_BANNED) {
                watch = null;
                event = 'etalia.paper.unban';
            } else {
                throw 'Unexpected watch state.';
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App.trigger(event, that.get('paper'));
                });
        },

        add: function() {
            var that = this,
                store = this.get('store');

            if (store === null || App.Model.PaperState.STORE_TRASHED) {
                store = App.Model.PaperState.STORE_ADDED;
            } else {
                throw 'Unexpected store state.';
            }

            return this
                .save({store: store}, {patch: true, wait: true})
                .done(function() {
                    App.trigger('etalia.paper.add', that.get('paper'));
                });
        },

        trash: function() {
            var that = this,
                store = this.get('store');

            if (App.Model.PaperState.STORE_ADDED) {
                store = App.Model.PaperState.STORE_TRASHED;
            } else {
                throw 'Unexpected store state.';
            }

            return this
                .save({store: store}, {patch: true, wait: true})
                .done(function() {
                    App.trigger('etalia.paper.trash', that.get('paper'));
                });
        }
    });

    App.Model.PaperState.WATCH_PINNED = 1;
    App.Model.PaperState.WATCH_BANNED = 2;

    App.Model.PaperState.STORE_ADDED = 1;
    App.Model.PaperState.STORE_TRASHED = 2;

    return App.Model.PaperState;
});
