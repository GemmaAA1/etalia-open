define([
    'app',
    'text!app/templates/thread/invite/treat.hbs',
    'app/model/thread/invite',
    'app/view/thread/thumb'
], function (App, template) {

    App.View.Thread.InviteTreat = App.Backbone.View.extend({

        tagName: 'div',
        className: 'invite-treat',

        template: App.Handlebars.compile(template),

        events: {
            "click .invite-treat-more": "onMoreClick",
            "click .invite-treat-accept": "onAcceptClick",
            "click .invite-treat-decline": "onDeclineClick"
        },

        initialize: function () {
            this.listenTo(this.model, "change", this.render);
        },

        onMoreClick: function(e) {
            e.preventDefault();

            var $button = this.$('.invite-treat-more'),
                $preview = this.$('.invite-treat-preview');
            if ($preview.hasClass('expanded')) {
                $preview.removeClass('expanded');
                $button.html('Read more&hellip;');
            } else {
                $preview.addClass('expanded');
                $button.html('Reduce');
            }
        },

        onAcceptClick: function(e) {
            e.preventDefault();

            var that = this;
            this.model
                .accept()
                .done(function() {
                    that.trigger('accepted');
                })
                .fail(function() {
                    App.log('Errors', 'Failed to accept invite');
                });
        },

        onDeclineClick: function(e) {
            e.preventDefault();

            var that = this;
            this.model
                .decline()
                .done(function() {
                    that.trigger('declined');
                })
                .fail(function() {
                    App.log('Errors', 'Failed to decline invite');
                });
        },

        render: function () {
            App.log('CommentThumbView::render');

            this.$el.html(this.template({
                content: this.model.get('thread').get('content')
            }));

            this.pushSubView(
                App.View.Thread.Thumb.create({
                    mode: App.View.Thread.Thumb.MODE_STANDALONE,
                    model: this.model.get('thread')
                }, {
                    $target: this.$('[data-thread-thumb]')
                })
            );

            return this;
        }
    });

    App.View.Thread.InviteTreat.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'options.model is expected to be an instance of App.Model.Invite';
        }
        if (!options.id) {
            options.id = 'thread-invite-treat-' + options.model.get('id');
        }

        var view = new App.View.Thread.InviteTreat(options);

        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Thread.InviteTreat;
});
