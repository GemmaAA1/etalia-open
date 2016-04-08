define([
    'app',
    'text!app/templates/user/thumb.html'
], function (App, template) {

    App.View.User = App.View.User || {};

    App.View.User.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'user',

        template: App.Handlebars.compile(template),

        // TODO click img : onImgClick (add 'hover' class)
        /*events: {
            "click .thumb-pin": "onPinClick",
            "click .thumb-ban": "onBanClick"
        },*/

        initialize: function () {
            this.listenTo(this.model, "change", this.render);
        },

        render: function() {
            this.$el.html(this.template(this.model.attributes));

            // TODO this.$el.addClass('followed');

            return this;
        }

        /*onPinClick: function(e) {
            e.preventDefault();
        },

        onBanClick: function(e) {
            e.preventDefault();
        }*/
    });

    App.View.User.Thumb.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'options.model is expected to be an instance of App.Model.User';
        }
        if (!options.id) {
            options.id = 'user-thumb-' + options.model.get('id');
        }

        var thumb = new App.View.User.Thumb(options);
        if (createOptions) {
            App.View.create(thumb, createOptions);
        }

        return thumb;
    };

    return App.View.User.Thumb;
});
