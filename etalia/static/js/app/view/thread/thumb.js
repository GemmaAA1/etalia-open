define([
    'app/app',
    'text!app/templates/thread/thumb.hbs',
    'app/view/user/thumb'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thumb',

        template: App.Handlebars.compile(template),
        buttons: {
            pin: false,
            ban: false,
            join: false,
            leave: false
        },

        mode: null,
        list: null,

        events: {
            "click .title a": "onTitleClick",
            "click .thumb-pin": "onPinClick",
            "click .thumb-ban": "onBanClick",
            "click .thumb-join": "onJoinClick",
            "click .thumb-leave": "onLeaveClick"
        },

        initialize: function (options) {
            this.mode = options.mode || App.View.Thread.Thumb.MODE_LIST;
            if (this.mode == App.View.Thread.Thumb.MODE_LIST) {
                if (!options.list) {
                    throw '"options.list" is mandatory is list mode.';
                }
                this.list = options.list;
            }

            if (options.buttons) {
                this.buttons = App._.extend({}, this.buttons, options.buttons);
            }

            this.listenTo(this.model, "change", this.render);
            this.listenTo(this.model, "change:state", this.render);
            this.listenTo(this.model, "add:posts remove:posts", this.updatePostsCount);
            this.listenTo(this.model, "add:members remove:members", this.updateMembersCount);
        },

        onTitleClick: function(e) {
            e.preventDefault();

            if (this.mode == App.View.Thread.Thumb.MODE_LIST) {
                this.list.trigger('model:detail', this.model, this);
            }
        },

        onPinClick: function(e) {
            e.preventDefault();

            this.model.getState().togglePinned();
        },

        onBanClick: function(e) {
            e.preventDefault();

            this.model.getState().toggleBanned();
        },

        onJoinClick: function(e) {
            e.preventDefault();

            this.model.getState().join();
        },

        onLeaveClick: function(e) {
            e.preventDefault();

            this.model.getState().leave();
        },

        updateMembersCount: function() {
            this.$('.icons .member .count').text(this.model.getMembersCount());
        },

        updatePostsCount: function() {
            this.$('.icons .comment .count').text(this.model.getPostsCount());
        },

        render: function() {
            App.log('ThreadThumbView::render');

            var //is_owner = this.model.isOwner(App.getCurrentUser()),
                is_member = this.model.isMember(App.getCurrentUser()),
                is_public = this.model.isPublic();

            var attributes = App._.extend({}, this.model.attributes, {
                list_mode: this.mode == App.View.Thread.Thumb.MODE_LIST,

                pin_button: this.buttons.pin,
                ban_button: this.buttons.ban,
                join_button: this.buttons.join && is_public && !is_member,
                leave_button: this.buttons.leave && is_member,

                members_count: this.model.getMembersCount(),
                posts_count: this.model.getPostsCount()
            });

            this.$el.html(this.template(attributes));

            this.pushSubView(
                App.View.User.Thumb.create({
                    model: this.model.get('user')
                }, {
                    $target: this.$('[data-user-placeholder]')
                })
            );

            return this;
        }
    });

    App.View.Thread.Thumb.MODE_LIST = 1;
    App.View.Thread.Thumb.MODE_STANDALONE = 2;

    App.View.Thread.Thumb.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Thread.Thumb(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Thread.Thumb;
});
