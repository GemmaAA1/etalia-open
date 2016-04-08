define([
    'app',
    'text!app/templates/thread/post/list.html',
    'app/view/thread/post/thumb',
    'app/view/thread/post/form'
], function (App, template) {

    App.View.Thread.PostList = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-posts',

        template: App.Handlebars.compile(template),

        thread: null,
        form: null,

        /*events: {
         "click .title a": "onTitleClick",
         "click .thumb-pin": "onPinClick",
         "click .thumb-ban": "onBanClick"
         },*/

        initialize: function (options) {
            if (!options.thread) {
                throw 'Expected instance of App.Model.Thread';
            }

            this.thread = options.thread;

            //this.listenTo(this.model, "change:state", this.onThreadStateChange);
            this.listenTo(this.model, "sync", this.render);
        },

        renderForm: function () {
            if (this.form) {
                this.form.remove();
            }

            this.form = App.View.Thread.PostForm.create({}, {
                $target: this.$('[data-post-add-form]')
            });
            /*this.$('.thread-post-form').html(
             this.form.render().$el
             );*/

            this.listenToOnce(this.form, 'validation_success', this.submitForm);
        },

        submitForm: function () {
            var that = this;
            this.form.model.save({
                user: App.Model.User.getCurrent(),
                thread: this.thread
            }, {
                success: function () {
                    that.renderForm();
                }
            });
        },

        render: function () {
            App.log('PostListView::render');

            this.$el.html(this.template()); // this.model.toJSON()

            var $list = this.$('.thread-posts-list');
            this.model.each(function (post) {
                App.View.Thread.PostThumb.create({
                    model: post
                }, {
                    $target: $list,
                    append: true
                });
            });

            App.View.User.Thumb.create({
                model: App.Model.User.getCurrent()
            }, {
                $target: this.$('.user-placeholder')
            });

            this.renderForm();

            return this;
        }
    });

    App.View.Thread.PostList.create = function (options, createOptions) {
        options = options || {};
        if (!options.thread) {
            throw 'options.thread is expected to be instance of App.Model.Thread';
        }
        options.model = options.thread.get('posts');

        var list = new App.View.Thread.PostList(options);

        if (createOptions) {
            App.View.create(list, createOptions);
        }

        return list;
    };

    return App.View.Thread.PostList;
});
