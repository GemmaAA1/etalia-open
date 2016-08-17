define(['app/app'], function (App) {

    var userPopoverPath = '/popover/states/';

    App.Model.Popover = App.Backbone.RelationalModel.extend({
        defaults: {
            title: null,
            body: null,
            anchor: null,
            type: null
        }
    });

    App.Model.Popover.TYPE_ANCHORED = 1;
    App.Model.Popover.TYPE_MODAL = 2;


    App.Model.UserPopover = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + userPopoverPath,

        defaults: {
            status: null,
            display: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'user',
                relatedModel: App.Model.User,
                includeInJSON: false,
                autoFetch: false
            },
            {
                type: App.Backbone.HasOne,
                key: 'popover',
                relatedModel: App.Model.Popover,
                includeInJSON: false,
                autoFetch: false
            }
        ],

        markDone: function() {
            var that = this,
                status = this.get('status');

            if (status === App.Model.UserPopover.STATUS_NEW) {
                status = App.Model.UserPopover.STATUS_DONE;
            } else {
                throw 'Unexpected user popover state.';
            }

            return this
                .save({status: status}, {patch: true, wait: true})
                .done(function() {
                    App.trigger('etalia.user_popover.mark_done', that);
                });
        }
    });

    App.Model.UserPopover.STATUS_NEW = 1;
    App.Model.UserPopover.STATUS_DONE = 2;
    App.Model.UserPopover.DISPLAY_HIDE = 1;
    App.Model.UserPopover.DISPLAY_SHOW = 2;

    App.Model.UserPopovers = App.Backbone.Collection.extend({
        url: App.config.api_root + userPopoverPath,
        model: App.Model.UserPopover,

        initialize: function(options) {
            this.setQuery(options ? options.query : {});
        },

        setQuery: function(query) {
            this.queryParams = App._.extend({
                status: App.Model.UserPopover.STATUS_NEW,
                display: App.Model.UserPopover.DISPLAY_SHOW
            }, query);

            return this;
        }
    });

    return App.Model.UserPopover;
});
