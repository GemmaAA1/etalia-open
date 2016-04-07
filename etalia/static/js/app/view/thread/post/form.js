define([
    'underscore',
    'app',
    'app/model/thread/post'
], function (_, App) {

    App.View.Thread.PostForm = App.Backbone.Form.extend({

        schema: {
            content: {type: 'Tinymce', validators: ['required'], title: false}
        },

        template: _.template('\
            <form class="thread-post-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <button type="submit" class="btn btn-primary thread-post-form-submit">\
                    Post\
                </button>\
            </form>\
        '),

        events: {
            "click .thread-post-form-submit": "onSubmit",
            "click .thread-post-form-cancel": "onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            if (0 < this.model.get('id')) {
                this.template = _.template('\
                    <form class="thread-post-form form-horizontal" role="form">\
                        <div data-fieldsets></div>\
                        <button type="submit" class="btn btn-primary thread-post-form-submit">\
                            OK\
                        </button>\
                        <button type="button" class="btn btn-default thread-post-form-cancel">\
                            Cancel\
                        </button>\
                    </form>\
                ');
            }

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

    App.View.Thread.PostForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = new App.Model.Post();
        }

        var form = new App.View.Thread.PostForm(options);

        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.PostForm;
});
