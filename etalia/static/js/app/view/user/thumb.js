define([
    'app/app',
    'text!app/templates/user/thumb.hbs'
], function (App, template) {

    App.View.User = App.View.User || {};

    App.View.User.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'user',

        template: App.Handlebars.compile(template),

        // TODO click img : onImgClick (add 'hover' class)
        events: {
            "click .user-follow": "onFollowClick",
            "click .user-unfollow": "onUnFollowClick"
        },

        initialize: function () {
            this.listenTo(this.model, "change", this.render);
        },

        render: function() {
            //App.log('UserThumbView::render', this.model.get('id'));

            var that = this;

            App.getCurrentUser()
                .isFollowed(that.model)
                .then(function (followed) {
                    var attributes = App._.extend({}, that.model.attributes, {
                        followed: followed,
                        current: App.getCurrentUser().get('id') === that.model.get('id')
                    });

                    that.$el.html(that.template(attributes));

                    if (followed) {
                        that.$el.addClass('followed');
                    } else {
                        that.$el.removeClass('followed');
                    }
                });

            return this;
        },

        onFollowClick: function(e) {
            e.preventDefault();

            App.getCurrentUser().follow(this.model);
        },

        onUnFollowClick: function(e) {
            e.preventDefault();

            App.getCurrentUser().unfollow(this.model);
        }
    });

    App.View.User.Thumb.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'options.model is expected to be an instance of App.Model.User';
        }
        if (!options.id) {
            options.id = 'user-thumb-' + options.model.get('id');
        }

        var view = new App.View.User.Thumb(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.User.Thumb;
});
