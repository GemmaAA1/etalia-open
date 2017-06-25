define([
    'app',
    'app/model/user/user',
    'app/model/user/user-lib'
], function (App) {

    var urlRoot = App.config.api_root + '/library/states';

    App.Model.PaperState = App.Backbone.RelationalModel.extend({

        urlRoot: urlRoot,

        defaults: {
            watch: null,
            store: null,
            is_orcid: null
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
                events = [];

            if (watch === App.Model.ThreadState.WATCH_BANNED) {
                events.push('etalia.paper.unban');
            }

            if (watch === App.Model.PaperState.WATCH_PINNED) {
                watch = null;
                events.push('etalia.paper.unpin');
            } else {
                watch = App.Model.PaperState.WATCH_PINNED;
                events.push('etalia.paper.pin');
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App._.each(events, function(event) {
                        App.trigger(event, that.get('paper'));
                    });
                });
        },

        toggleBanned: function() {
            var that = this,
                watch = this.get('watch'),
                events = [];

            if (watch === App.Model.ThreadState.WATCH_PINNED) {
                events.push('etalia.paper.unpin');
            }

            if (watch === App.Model.PaperState.WATCH_BANNED) {
                watch = null;
                events.push('etalia.paper.unban');
            } else {
                watch = App.Model.PaperState.WATCH_BANNED;
                events.push('etalia.paper.ban');
            }

            return this
                .save({watch: watch}, {patch: true, wait: true})
                .done(function() {
                    App._.each(events, function(event) {
                        App.trigger(event, that.get('paper'));
                    });
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

    App.Model.PaperState.emptyTrash = function () {
        return App.$.ajax({
            url: urlRoot + '/empty-trash',
            method: 'DELETE',
            dataType: 'json'
        });
    };

    App.Model.PaperState.WATCH_PINNED = 1;
    App.Model.PaperState.WATCH_BANNED = 2;

    App.Model.PaperState.STORE_ADDED = 1;
    App.Model.PaperState.STORE_TRASHED = 2;

    return App.Model.PaperState;
});
