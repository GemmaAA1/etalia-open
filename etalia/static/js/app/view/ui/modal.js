define([
    'app',
    'text!app/templates/ui/modal.hbs',
    'bootstrap'
], function (App, template) {

    var _ = App._,
        defaultButton = {
            label: 'Button',
            attr: {
                'type': 'button',
                'class': 'btn btn-default'
            }
        },
        closeButton = {
            label: 'Close',
            attr: {
                'type': 'button',
                'class': 'btn btn-default',
                'data-dismiss': 'modal'
            }
        };

    App.View.Ui = App.View.Ui || {};

    App.View.Ui.Modal = App.Backbone.View.extend({

        template: App.Handlebars.compile(template),

        buttonTemplate: App.Handlebars.compile(
            '<button{{#each attr}} {{@key}}="{{this}}"{{/each}}>{{label}}</button>'
        ),

        attributes: {
            'class': 'modal fade',
            'tabindex': '-1',
            'role': 'dialog'
        },

        defaults: {
            title: null,
            content: null,
            backdrop: true,
            footer: {
                buttons: [closeButton]
            }
        },

        events: {
            "hidden.bs.modal": "_onHidden",
            "hide.bs.modal": "_onHide"
        },

        initialize: function (options) {
            options || (options = {});
            _.defaults(this, this.defaults);
            _.extend(this, _.pick(options, _.keys(this.defaults)));

            _.bindAll(this, 'close', '_onHidden', '_onHide');
        },

        render: function () {
            var that = this,
                renderFooter = _.isObject(this.footer) && _.isObject(this.footer.buttons);

            this.$el.html(this.template({
                title: this.title,
                footer: renderFooter
            }));

            // Renders the content
            if (_.isObject(this.content) && _.isFunction(this.content.render)) {
                this.$('.modal-body').append(this.content.render().$el);
            } else {
                this.$('.modal-body').append(this.content);
            }


            // Renders the footer
            if (renderFooter) {
                var $footer = this.$('.modal-footer');
                _.each(this.footer.buttons, function (button) {
                    _.defaults(button, defaultButton);
                    $footer.append(that.buttonTemplate(button));
                });
            }

            // Call bootstrap modal
            this.$el.modal({
                keyboard: false,
                backdrop: this.backdrop
            });
        },

        _onHidden: function () {
            if (_.isObject(this.content) && _.isFunction(this.content.remove)) {
                this.content.remove();
            }
            this.remove();

            this.trigger('hidden');
        },

        _onHide: function () {
            this.trigger('hide');
        },

        close: function () {
            this.$el.modal("hide");

            return this;
        }
    });

    App.View.Ui.Modal.defaultButton = closeButton;

    return App.View.Ui.Modal;
});
