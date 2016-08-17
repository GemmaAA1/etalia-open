define([
    'app/app',
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
        },
        defaults = {
            title: null,
            content: null,
            backdrop: true,
            footer: {
                buttons: [closeButton]
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
            //'tabindex': '-1', Select2 bug : http://stackoverflow.com/questions/18487056/select2-doesnt-work-when-embedded-in-a-bootstrap-modal
            'role': 'dialog'
        },

        events: {
            "shown.bs.modal": "_onShown",
            "show.bs.modal": "_onShow",
            "hidden.bs.modal": "_onHidden",
            "hide.bs.modal": "_onHide",
            "click .modal-footer button": "_onFooterButtonClick"
        },

        initialize: function (options) {
            options || (options = {});
            _.defaults(this, defaults);
            _.extend(this, _.pick(options, _.keys(defaults)));

            _.bindAll(this, 'close', '_onHidden', '_onHide');
        },

        render: function () {
            this.$el.html(this.template({
                title: this.title
            }));

            this._renderContent();
            this._renderFooter();

            // Call bootstrap modal
            this.$el.modal({
                keyboard: false,
                backdrop: this.backdrop
            });

            return this;
        },

        updateTitle: function(title) {
            this.$('.modal-title').html(title);

            return this;
        },

        updateContent: function (content) {
            if (!content) {
                throw 'content argument is mandatory';
            }
            this._removeContent();
            this.content = content;
            this._renderContent();

            return this;
        },

        updateFooter: function(footer) {
            this.footer = footer;
            this._renderFooter();

            return this;
        },

        _removeContent: function() {
            if (_.isObject(this.content) && _.isFunction(this.content.remove)) {
                this.content.remove();
            } else {
                this.$('.modal-body').empty();
            }
        },

        _renderContent: function() {
            if (_.isObject(this.content) && _.isFunction(this.content.render)) {
                this.$('.modal-body').append(this.content.render().$el);
            } else {
                this.$('.modal-body').append(this.content);
            }
        },

        _renderFooter: function() {
            var that = this,
                $footer = this.$('.modal-footer').empty();
            if (_.isObject(this.footer) && _.isObject(this.footer.buttons)) {
                $footer.show();
                _.each(this.footer.buttons, function (button) {
                    _.defaults(button, defaultButton);
                    $footer.append(that.buttonTemplate(button));
                });
            } else {
                $footer.hide();
            }
        },

        _onShown: function () {
            this.trigger('shown');
        },

        _onShow: function () {
            this.trigger('show');
        },

        _onHidden: function () {
            this._removeContent();
            this.remove();

            this.trigger('hidden');
        },

        _onHide: function () {
            this.trigger('hide');
        },

        _onFooterButtonClick: function(e) {
            this.trigger('button_click', e);
        },

        close: function () {
            this.$el.modal("hide");

            return this;
        }
    });

    App.View.Ui.Modal.defaultButton = closeButton;

    return App.View.Ui.Modal;
});
