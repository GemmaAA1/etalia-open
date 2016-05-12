define([
    'app',
    'text!app/templates/thread/thumb.hbs',
    'app/view/user/thumb'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    App.View.Thread.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thumb',

        template: App.Handlebars.compile(template),

        mode: null,
        list: null,

        events: {
            "click .title a": "onTitleClick",
            "click .thumb-pin": "onPinClick",
            "click .thumb-ban": "onBanClick"
        },

        initialize: function (options) {
            this.mode = options.mode || App.View.Thread.Thumb.MODE_LIST;
            if (this.mode == App.View.Thread.Thumb.MODE_LIST) {
                if (!options.list) {
                    throw '"options.list" is mandatory is list mode.';
                } else {
                    this.list = options.list;
                }
            }

            this.listenTo(this.model, "change", this.render);
            this.listenTo(this.model, "change:state", this.onThreadStateChange);
            this.listenTo(this.model, "add:posts remove:posts", this.updatePostsCount);
            this.listenTo(this.model, "add:members remove:members", this.updateMembersCount);
        },

        onTitleClick: function(e) {
            e.preventDefault();

            if (this.mode == App.View.Thread.Thumb.MODE_LIST) {
                this.list.trigger('model:detail', this.model, this);
            }
        },

        onThreadStateChange: function(state, thread) {
            if (state.get('watch') === App.Model.State.WATCH_BANNED) {
                this.list.trigger('model:remove', thread);
            } else {
                this.render();
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

        updateMembersCount: function() {
            this.$('.icons .member .count').text(this.model.getMembersCount());
        },

        updatePostsCount: function() {
            this.$('.icons .comment .count').text(this.model.getPostsCount());
        },

        render: function() {
            App.log('ThreadThumbView::render');

            var attributes = App._.extend(this.model.attributes, {
                list_mode: this.mode == App.View.Thread.Thumb.MODE_LIST,
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
