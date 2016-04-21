define([
    'underscore',
    'app',
    'app/model/thread/post'
], function (_, App) {

    App.View.Thread.EditForm = App.Backbone.Form.extend({

        schema: {
            content: {type: 'Tinymce', title: false, settings: {height: 400}}
        },

        template: _.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-footer">\
                    <button type="submit" class="btn btn-primary form-submit">\
                        Update\
                    </button>\
                    <button type="button" class="btn btn-default form-cancel">\
                        Cancel\
                    </button>\
                </div>\
            </form>\
        '),

        events: {
            "click .form-submit": "onSubmit",
            "click .form-cancel": "onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            this.listenTo(this, "submit", this.onSubmit);
        },

        onSubmit: function(e) {
            e.preventDefault();

            var errors = this.commit({ validate: true });
            if (errors) {
                App.log('Errors', errors);

                this.trigger('validation_failure', this);
                return;
            }

            this.trigger('validation_success', this);
        },

        onCancel: function(e) {
            e.preventDefault();

            this.trigger('cancel', this);
        }
    });


    App.View.Thread.EditForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'Expected options.model instance of App.Model.Thread';
        }

        var form = new App.View.Thread.EditForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.EditForm;
});
