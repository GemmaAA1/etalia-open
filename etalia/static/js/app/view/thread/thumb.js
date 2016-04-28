define([
    'app',
    'text!app/templates/thread/thumb.html',
    'app/view/user/thumb'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thumb',

        template: App.Handlebars.compile(template),

        list: null,

        events: {
            "click .title a": "onTitleClick",
            "click .thumb-pin": "onPinClick",
            "click .thumb-ban": "onBanClick"
        },

        initialize: function (options) {
            if (!options.list) {
                throw '"options.list" is mandatory.';
            }
            this.list = options.list;

            this.listenTo(this.model, "change:state", this.onThreadStateChange);
            this.listenTo(this.model, "change", this.render);
            this.listenTo(this.model, "add:posts remove:posts", this.updatePostsCount);
            this.listenTo(this.model, "add:members remove:members", this.updateMembersCount);
        },

        onTitleClick: function(e) {
            e.preventDefault();

            this.list.trigger('model:detail', this.model, this);
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

            this.model.getState()
                .togglePinned()
                .save(null, {wait:true});
        },

        onBanClick: function(e) {
            e.preventDefault();

            var thread = this.model;

            thread.getState()
                .toggleBanned()
                .save(null, {wait:true});
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
});
