define([
    'app',
    'text!app/templates/paper/thumb.hbs'
], function (App, template) {

    App.View.Paper = App.View.Paper || {};

    App.View.Paper.Thumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thumb',

        template: App.Handlebars.compile(template),
        buttons: {
            pin: false,
            ban: false,
            add: false,
            trash: false
        },

        mode: null,
        list: null,

        events: {
            "click .title a": "onTitleClick",
            "click .thumb-pin": "onPinClick",
            "click .thumb-ban": "onBanClick",
            "click .thumb-add": "onAddClick",
            "click .thumb-trash": "onTrashClick"
        },

        initialize: function (options) {
            this.mode = options.mode || App.View.Paper.Thumb.MODE_LIST;
            if (this.mode == App.View.Paper.Thumb.MODE_LIST) {
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
        },

        onTitleClick: function(e) {
            e.preventDefault();

            if (this.mode == App.View.Paper.Thumb.MODE_LIST) {
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

        onAddClick: function(e) {
            e.preventDefault();

            this.model.getState().add();
        },

        onTrashClick: function(e) {
            e.preventDefault();

            this.model.getState().trash();
        },

        render: function() {
            //App.log('PaperThumbView::render');

            var attributes = App._.extend({}, this.model.attributes, {
                list_mode: this.mode == App.View.Paper.Thumb.MODE_LIST,

                pin_button: this.buttons.pin,
                ban_button: this.buttons.ban,
                add_button: this.buttons.add && !this.model.isAdded(),
                trash_button: this.buttons.trash && this.model.isAdded()
            });

            this.$el.html(this.template(attributes));

            return this;
        }
    });

    App.View.Paper.Thumb.MODE_LIST = 1;
    App.View.Paper.Thumb.MODE_STANDALONE = 2;

    App.View.Paper.Thumb.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Paper.Thumb(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Paper.Thumb;
});
