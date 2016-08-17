define(['app/app', 'app/model/thread/comment'], function (App) {

    App.View.Thread.CommentForm = App.Backbone.Form.extend({

        schema: {
            content: {type: 'TextArea', validators: ['required'], title: 'Your comment'}
        },

        template: _.template('\
            <form class="thread-comment-form form-horizontal" role="form">\
                <button type="submit" class="btn btn-primary thread-comment-form-submit">\
                    <span class="eai eai-tick"></span>\
                </button>\
                <button type="button" class="btn btn-primary thread-comment-form-cancel">\
                    <span class="eai eai-remove"></span>\
                </button>\
                <div data-fieldsets></div>\
            </form>\
        '),

        events: {
            "click .thread-comment-form-submit": "onSubmit",
            "click .thread-comment-form-cancel": "onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            this.listenTo(this, "submit", this.onSubmit);
        },

        onSubmit: function(e) {
            e.preventDefault();

            App.log('CommentFormView::onSubmit');

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

    App.View.Thread.CommentForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = new App.Model.Comment();
        }

        var form = new App.View.Thread.CommentForm(options);

        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.CommentForm;
});
