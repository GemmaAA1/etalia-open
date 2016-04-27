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
            this.listenTo(this.model, "sync change", this.render);
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
                .save();
        },

        onBanClick: function(e) {
            e.preventDefault();

            var thread = this.model;

            thread.getState()
                .toggleBanned()
                .save()
                .then(function() {
                    if (thread.get('state').get('watch') == App.Model.State.WATCH_BANNED) {

                    }
                });
        },

        render: function() {
            App.log('ThreadThumbView::render');

            this.$el.html(this.template(this.model.attributes));

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
